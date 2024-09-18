# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations

import typing as ty
from contextlib import contextmanager

from attrs import frozen
from loguru import logger
from sqlitedict import SqliteDict


@frozen
class CachedMapping:
    """A `CachedMapping` is a way to store a large number of objects in an
    sqlite database. This is basically a wrapper around an `SqliteDict` object
    with a few connivances for the typical usage patterns.
    """

    filename: str
    commit_size: int
    tables: ty.Set[str]

    @classmethod
    def empty(cls, filename="cached.sqlite", commit_size=100000) -> CachedMapping:
        """Create an empty `CachedMapping`.

        :commit_size: Number of objects to commit in one batch. This should
        generally be some large value as backing sqlite score can work with
        large values and small ones lead to serious performance issues
        :filename: Name of the sqlite file to use.
        """
        return cls(filename=filename, commit_size=commit_size, tables=set())

    def cache(self, name: str, entries: ty.Iterable[ty.Tuple[str, ty.Any]]):
        """Cache the values in the given iterable for later use.

        The iterable should provide a tuple of (key, value) pairs and each key
        should be a string and each value a Pickleable object.

        :name: Name of the cache.
        :results: The iterable to cache.
        """

        added = 0
        logger.debug("Caching values in {}.{}", self.filename, name)
        with SqliteDict(
            self.filename, tablename=name, autocommit=False, flag="c"
        ) as db:
            for key, value in entries:
                db[key] = value
                added += 1
                if added % self.commit_size == 0:
                    logger.trace("Committing to {}.{}", self.filename, name)
                    db.commit()
            logger.trace("Final commit to {}.{}", self.filename, name)
            db.commit()
            db.close()
        if added:
            self.tables.add(name)
            logger.debug("Cached {} values", added)
        else:
            logger.info("No values cached for {}", name)

    @contextmanager
    def cached_values(self, name: str, **kwargs) -> ty.Iterator[ty.Dict[str, ty.Any]]:
        """Get a context handler with cached values for the given cache.

        :name: The name of the cache to use.
        """
        if name not in self.tables:
            yield {}
        else:
            with SqliteDict(self.filename, tablename=name, flag="r", **kwargs) as db:
                yield db
