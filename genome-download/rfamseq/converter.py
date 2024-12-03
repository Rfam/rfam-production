# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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


from attrs import fields, has
from cattrs import Converter
from cattrs.gen import make_dict_structure_fn, make_dict_unstructure_fn, override


def to_camel_case(snake_str: str) -> str:
    components = snake_str.split("_")
    return components[0] + "".join(x.title() for x in components[1:])


def camel_case_converter() -> Converter:
    """A catters converter that properly handles snake case and camel case
    properties."""
    converter = Converter()

    def to_camel_case_unstructure(cls):
        return make_dict_unstructure_fn(
            cls,
            converter,
            **{a.name: override(rename=to_camel_case(a.name)) for a in fields(cls)}
        )

    def to_camel_case_structure(cls):
        return make_dict_structure_fn(
            cls,
            converter,
            **{a.name: override(rename=to_camel_case(a.name)) for a in fields(cls)}
        )

    converter.register_unstructure_hook_factory(has, to_camel_case_unstructure)
    converter.register_structure_hook_factory(has, to_camel_case_structure)
    return converter
