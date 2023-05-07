# -*- coding: utf-8 -*-

"""
Copyright [2009-2023] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import logging
import os
import subprocess as sp
import tempfile
import typing as ty
from contextlib import contextmanager
from pathlib import Path

from rfamseq.accession import Accession
from rfamseq.utils import assert_never
from rfamseq.wgs import WgsSequenceId

LOGGER = logging.getLogger(__name__)


@contextmanager
def id_handle(ids: ty.Iterable[str]) -> ty.Iterator[ty.IO]:
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8") as tmp:
        for id in ids:
            tmp.write(id)
            tmp.write("\n")
        tmp.seek(0)
        yield tmp


def index(path: Path):
    LOGGER.info("Indexing %s", path)
    sp.check_call(
        ["esl-sfetch", "--index", str(path)], stderr=sp.DEVNULL, stdout=sp.DEVNULL
    )


@contextmanager
def filtered(path: Path, ids: ty.Iterable[str]) -> ty.Iterator[ty.IO]:
    indexed = path.parent / (path.name + ".ssi")
    if not indexed.exists():
        index(path)
    with id_handle(ids) as to_select:
        with tempfile.NamedTemporaryFile(dir=os.curdir) as tmp_fasta:
            LOGGER.debug("Running esl-sfetch")
            sp.check_call(
                ["esl-sfetch", "-o", tmp_fasta.name, "-f", str(path), to_select.name],
                stderr=sp.DEVNULL,
                stdout=sp.DEVNULL,
            )
            LOGGER.debug("esl-sfetch finished")
            yield tmp_fasta
