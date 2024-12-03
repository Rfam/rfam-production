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


class UnknownFileType(Exception):
    """Raised when a file with an unhandleable filetype is downloaded"""


@contextmanager
def wget(url: str) -> ty.Iterator[ty.IO]:
    """Fetch the data at the given URL. This will make a reques tto the site
    and write the results to a temporary file that can be read. The file will
    be stored on disk to try and limit the memory usage. Additionally, the data
    will be decompressed if needed.
    """

    if not url:
        raise ValueError("Must give a url like string")

    # This uses the current directory because on the cluster the temp
    # space may be highly limited. We already run things on relatively
    # fast storage so it doesn't slow things down too much to use the
    # current directory for temp space.
    with tempfile.NamedTemporaryFile("wb+", dir=os.curdir) as tmp:
        LOGGER.debug("Fetching %s to %s", url, tmp.name)
        try:
            with closing(urllib.request.urlopen(str(url))) as req:
                shutil.copyfileobj(req, tmp, length=16 * 1024 * 1024)
        except Exception as err:
            LOGGER.exception(err)
            LOGGER.warning("Failed fetching %s", url)
            raise err
        tmp.flush()
        tmp.seek(0)

        info = sp.check_output(["file", "--mime-type", tmp.name])
        if b"application/gzip" in info:
            with tempfile.NamedTemporaryFile("w+", dir=os.curdir) as decomp:
                LOGGER.debug("Decompressing file for %s to %s", url, decomp.name)
                try:
                    sp.run(
                        ["zcat", "-f", tmp.name],
                        check=True,
                        stdout=decomp,
                        stderr=sp.PIPE,
                    )
                except sp.CalledProcessError as err:
                    LOGGER.warn("Failed to decrompress file %s", tmp.name)
                    LOGGER.warn(err.stderr)
                    LOGGER.exception(err)
                    raise err
                decomp.flush()
                decomp.seek(0)
                yield decomp
        elif b"text/plain" in info:
            with io.open(tmp.name, "r") as raw:
                yield raw
        else:
            raise UnknownFileType(f"Could not handle file type from {url} ({info})")
