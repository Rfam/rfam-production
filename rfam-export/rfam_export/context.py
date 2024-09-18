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

import typing as ty
from contextlib import contextmanager

import pymysql.cursors
from attr import Factory, frozen

from rfam_export.cached_mapping import CachedMapping


@frozen
class Context:
    rfam_version: str
    rfamseq_version: str
    mysql_host: str
    mysql_user: str
    mysql_password: str
    mysql_port: int
    mysql_database: str
    cached_mapping: CachedMapping = Factory(CachedMapping.empty)

    @contextmanager
    def connection(self):
        connection = pymysql.connect(
            host=self.mysql_host,
            user=self.mysql_user,
            password=self.mysql_password,
            port=self.mysql_port,
            database=self.mysql_database,
            cursorclass=pymysql.cursors.DictCursor,
        )
        yield connection
        connection.close()

    @contextmanager
    def cursor(self):
        with self.connection() as conn:
            with conn.cursor() as cursor:
                yield cursor

    def cache(self, name: str, entries: ty.Iterable[ty.Tuple[str, ty.Any]]):
        """Cache the given iterable of (key, value) pairs in an sqlite database.

        :name: Name of the cache.
        :entries: An Iterable of (key, value) pairs to store.
        """
        self.cached_mapping.cache(name, entries)

    @contextmanager
    def cached(self, name: str, **kwargs) -> ty.Iterator[ty.Dict[str, ty.Any]]:
        with self.cached_mapping.cached_values(name, **kwargs) as cache:
            yield cache
