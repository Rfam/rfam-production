# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
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
import shutil
import subprocess as sp
import tempfile
import typing as ty
import urllib.request
from contextlib import closing, contextmanager

LOGGER = logging.getLogger(__name__)


@contextmanager
def wget(url: str) -> ty.Iterator[ty.IO]:
    with tempfile.NamedTemporaryFile("wb+", dir=os.curdir) as tmp:
        LOGGER.debug("Fetching %s to %s", url, tmp.name)
        try:
            with closing(urllib.request.urlopen(str(url))) as req:
                shutil.copyfileobj(req, tmp, length=16 * 1024 * 1024)
        except Exception as err:
            LOGGER.warn("Failed fetching %s", url)
            raise err
        tmp.flush()
        tmp.seek(0)

        info = sp.check_output(["file", "--mime-type", tmp.name])
        if b"application/gzip" in info:
            with tempfile.NamedTemporaryFile("w+", dir=os.curdir) as decomp:
                LOGGER.info("Decompressing file for %s to %s", url, decomp.name)
                try:
                    sp.run(["zcat", "-f", tmp.name], check=True, stdout=decomp)
                except sp.CalledProcessError as err:
                    LOGGER.warn("Failed to decrompress file %s", tmp.name)
                    raise err
                decomp.flush()
                decomp.seek(0)
                yield decomp
        elif b"text/plain" in info:
            with io.open(tmp.name, "r") as raw:
                yield raw
        else:
            raise ValueError(f"Could not handle file type from {url} ({info})")
