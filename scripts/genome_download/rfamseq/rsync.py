# -*- coding: utf-8 -*-

"""
Copyright [2009-${2023}] EMBL-European Bioinformatics Institute
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


import io
import logging
import os
import subprocess as sp
import tempfile
import typing as ty
from contextlib import contextmanager
from pathlib import Path

LOGGER = logging.getLogger(__name__)


@contextmanager
def rsync(path: Path) -> ty.Iterator[ty.IO]:
    """
    Try to rsync a path to the local directory to read it. This will uncompress
    the file if needed.
    """

    with tempfile.NamedTemporaryFile("wb+", dir=os.curdir) as tmp:
        LOGGER.debug("Copying %s to %s", path, tmp.name)
        sp.check_call(["rsync", "-av", str(path), tmp.name])
        tmp.flush()

        info = sp.check_output(["file", "--mime-type", tmp.name])
        if b"application/gzip" in info:
            with tempfile.NamedTemporaryFile("w+", dir=os.curdir) as decomp:
                LOGGER.info("Decompressing file for %s to %s", path, decomp.name)
                sp.run(["zcat", "-f", tmp.name], check=True, stdout=decomp)
                decomp.flush()
                decomp.seek(0)
                yield decomp
        elif b"text/plain" in info:
            with io.open(tmp.name, "r") as raw:
                yield raw
        else:
            raise ValueError(f"Could not handle file type from {path} ({info})")
