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

from datetime import date, datetime
from itertools import islice
from typing import Iterable, List, NoReturn, TypeVar

T = TypeVar("T")


def assert_never(x: NoReturn) -> NoReturn:
    assert False, "Unhandled type: {}".format(type(x).__name__)


def serialize(obj):
    """
    Serialize datetime and date objects that are not serializable by default json code
    """

    if isinstance(obj, (datetime, date)):
        return obj.isoformat()
    raise TypeError("Type {obj} is not JSON serializable".format(obj=type(obj)))


def batched(iterable: Iterable[T], n: int) -> Iterable[List[T]]:
    "Batch data into tuples of length n. The last batch may be shorter."
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := list(islice(it, n)):
        yield batch
