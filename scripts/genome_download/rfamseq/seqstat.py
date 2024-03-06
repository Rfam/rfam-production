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

import re
import typing as ty

from attr import define, frozen


@frozen
class SeqStatInfo:
    """This represents the a sequence information entry from seqstat.

    :param sequence_id str: The sequence id of the entry.
    :param length int: The length of the sequence.
    :param description str: The description of the sequence.
    """

    sequence_id: str
    length: int
    description: str

    @classmethod
    def from_line(cls, line: str) -> SeqStatInfo:
        """Create a SeqStatInfo entry from a seqstat sequence line. These lines
        must start with an '='.

        >>> line = "= JADOXO010001145.1             2462 Postia placenta strain FPRL280 NODE_1146, whole genome shotgun sequence"
        >>> SeqStatInfo.from_line(line)
        SeqStatInfo(sequence_id="JADOXO010001145.1", length=2462, description="Postia placenta strain FPRL280 NODE_1146, whole genome shotgun sequence")
        """

        assert line[0] == "=", "Invalid seqstat line"
        parts = re.split(r"\s+", line, 3)
        return cls(sequence_id=parts[1], length=int(parts[2]), description=parts[3])

    @property
    def size_estimate(self) -> int:
        """This gives an estimate of the size of this sequence in bytes. This
        assumes that the sequence is ASCII characters and thus the byte size is
        the length, further is ignores the newlines the header (>id
        description) that would be included in a fasta file. This basically
        assumes that that sequence will be the largest, or most important, part of
        the written file.
        """
        return self.length


@define
class SeqStatCollection:
    """This is a collection of SeqStatInfo objects, that also tracks the
    current size. This should onl ybe modified with the `add` method.
    """

    info: ty.List[SeqStatInfo]
    current_size: float

    @classmethod
    def empty(cls) -> SeqStatCollection:
        """Create a new empty SeqStatCollection."""
        return SeqStatCollection(info=[], current_size=0.0)

    @property
    def entries(self) -> ty.List[SeqStatInfo]:
        """All entries in this collection."""
        return self.info

    def add(self, entry: SeqStatInfo):
        """Add a SeqStatInfo entry to this collection."""
        self.current_size += entry.size_estimate
        self.info.append(entry)

    @property
    def size_estimate(self) -> float:
        """The current size estimate of all entries."""
        return self.current_size


def parse(handle: ty.IO) -> ty.Iterator[SeqStatInfo]:
    """Parse the given file handle and produce an iterable of all SeqStatInfo
    entries in the file.
    """

    for line in handle:
        if not line.startswith("="):
            continue
        yield SeqStatInfo.from_line(line)


def chunk(
    iterable: ty.Iterator[SeqStatInfo], byte_size: int
) -> ty.Iterator[SeqStatCollection]:
    """This will chunk the given iterable of SeqStatInfo entries into chunks
    where each one is roughly byte_size. Size must be in bytes. The produced
    size will often be larger than the target size as it only breaks the group
    once it is larger than the given size.

    :param iterable: An iterable of SeqStatInfo entries.
    :param byte_size: The target size in bytes.
    """

    chunk = SeqStatCollection.empty()
    for info in iterable:
        chunk.add(info)
        if chunk.size_estimate > byte_size:
            yield chunk
            chunk = chunk.empty()
    yield chunk
