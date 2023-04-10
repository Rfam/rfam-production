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


import pytest

from rfamseq import wgs


@pytest.mark.parametrize(
    "kind,raw,expected",
    [
        (
            wgs.WgsSequenceKind.SCAFFOLD,
            "ALWZ04S0000001-ALWZ04S3033285",
            wgs.WgsSequence(
                wgs_id="ALWZ04S",
                kind=wgs.WgsSequenceKind.SCAFFOLD,
                start=1,
                stop=3033285,
                id_length=len("ALWZ04S0000001"),
            ),
        ),
        (
            wgs.WgsSequenceKind.SEQUENCE,
            "ALWZ040000001-ALWZ046221640",
            wgs.WgsSequence(
                wgs_id="ALWZ04",
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=1,
                stop=6221640,
                id_length=len("ALWZ040000001"),
            ),
        ),
        (
            wgs.WgsSequenceKind.SEQUENCE,
            "AMCR01000001-AMCR01022363",
            wgs.WgsSequence(
                wgs_id="AMCR01",
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=1,
                stop=22363,
                id_length=len("AMCR01000001"),
            ),
        ),
        (
            wgs.WgsSequenceKind.CONTIG,
            "ACTP02000001-ACTP02000023",
            wgs.WgsSequence(
                wgs_id="ACTP02",
                kind=wgs.WgsSequenceKind.CONTIG,
                start=1,
                stop=23,
                id_length=len("ACTP02000001"),
            ),
        ),
    ],
)
def test_wgs_info_parses_correctly(kind, raw, expected):
    info = wgs.WgsSequence.build(kind, raw)
    assert info == expected


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("KE150405-KE150411", wgs.ContigInfo(prefix="KE", start=150405, stop=150411)),
    ],
)
def test_can_parse_scaffolds(raw, expected):
    info = wgs.ContigInfo.build(raw)
    assert info == expected


@pytest.mark.parametrize(
    "kind,raw",
    [
        (wgs.WgsSequenceKind.SCAFFOLD, "ALWZ04S0000001-ALWZ04S3033285"),
        (wgs.WgsSequenceKind.SEQUENCE, "ALWZ040000001-ALWZ046221640"),
        (wgs.WgsSequenceKind.SEQUENCE, "AMCR01000001-AMCR01022363"),
    ],
)
def test_parsing_preseves_ids(kind, raw):
    info = wgs.WgsSequence.build(kind, raw)
    assert info.as_range() == raw


@pytest.mark.parametrize(
    "info,id,expected",
    [
        (
            wgs.WgsSequence(
                wgs_id="ALWZ04S",
                kind=wgs.WgsSequenceKind.SCAFFOLD,
                start=1,
                stop=3033285,
                id_length=len("ALWZ04S0000001"),
            ),
            "ALWZ04S0000001",
            True,
        ),
        (
            wgs.WgsSequence(
                wgs_id="ALWZ04S",
                kind=wgs.WgsSequenceKind.SCAFFOLD,
                start=1,
                stop=3033285,
                id_length=len("ALWZ04S0000001"),
            ),
            "ALWZ04S9999999",
            False,
        ),
        (
            wgs.WgsSequence(
                wgs_id="AMCR01",
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=1,
                stop=22363,
                id_length=len("AMCR01000001"),
            ),
            "AMCR01000001",
            True,
        ),
    ],
)
def test_wgs_knows_what_it_contains(info, id, expected):
    assert (id in info) == expected


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("ALWZ04S0000001", ("ALWZ04S", 1)),
        ("ALWZ04S6221640", ("ALWZ04S", 6221640)),
        ("ALWZ046221640", ("ALWZ04", 6221640)),
        ("LKAM01000000", ("LKAM01", 0)),
        ("AMCR01000000", ("AMCR01", 0)),
    ],
)
def test_can_get_endpoint_correctly(raw, expected):
    assert wgs.wgs_endpoint(raw) == expected


@pytest.mark.parametrize(
    "info,expected",
    [
        (
            wgs.WgsSequence(
                wgs_id="ALWZ04",
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=1,
                stop=10,
                id_length=len("ALWZ040000001"),
            ),
            [
                "ALWZ040000001",
                "ALWZ040000002",
                "ALWZ040000003",
                "ALWZ040000004",
                "ALWZ040000005",
                "ALWZ040000006",
                "ALWZ040000007",
                "ALWZ040000008",
                "ALWZ040000009",
                "ALWZ040000010",
            ],
        ),
        (
            wgs.WgsSequence(
                wgs_id="ALWZ04S",
                kind=wgs.WgsSequenceKind.SCAFFOLD,
                start=1,
                stop=10,
                id_length=len("ALWZ04S0000001"),
            ),
            [
                "ALWZ04S0000001",
                "ALWZ04S0000002",
                "ALWZ04S0000003",
                "ALWZ04S0000004",
                "ALWZ04S0000005",
                "ALWZ04S0000006",
                "ALWZ04S0000007",
                "ALWZ04S0000008",
                "ALWZ04S0000009",
                "ALWZ04S0000010",
            ],
        ),
        (
            wgs.WgsSequence(
                wgs_id="AMCR01",
                kind=wgs.WgsSequenceKind.SCAFFOLD,
                start=1,
                stop=10,
                id_length=len("AMCR010000001"),
            ),
            [
                "AMCR010000001",
                "AMCR010000002",
                "AMCR010000003",
                "AMCR010000004",
                "AMCR010000005",
                "AMCR010000006",
                "AMCR010000007",
                "AMCR010000008",
                "AMCR010000009",
                "AMCR010000010",
            ],
        ),
    ],
)
def test_wgs_generates_expected_ids(info, expected):
    assert list(info.ids()) == expected


@pytest.mark.parametrize(
    "accession,expected",
    [
        (
            "ALWZ05",
            (
                [],
                [
                    wgs.WgsSequence(
                        wgs_id="ALWZ05",
                        kind=wgs.WgsSequenceKind.SEQUENCE,
                        start=1,
                        stop=5154716,
                        id_length=13,
                    ),
                    wgs.WgsSequence(
                        wgs_id="ALWZ05S",
                        kind=wgs.WgsSequenceKind.CONTIG,
                        start=1,
                        stop=2066983,
                        id_length=len("ALWZ05S0000001"),
                    ),
                ],
            ),
        ),
        (
            "ACTP02",
            (
                [wgs.ContigInfo(prefix="KE", start=150405, stop=150411)],
                [
                    wgs.WgsSequence(
                        wgs_id="ACTP02",
                        kind=wgs.WgsSequenceKind.SEQUENCE,
                        start=1,
                        stop=23,
                        id_length=len("ACTP02000001"),
                    ),
                ],
            ),
        ),
        # (
        #     "SCMI01",
        #     [
        #         wgs.WgsSequence(
        #             wgs_id='SCMI01',
        #             kind=wgs.WgsSequenceKind.SEQUENCE,
        #             start=1,
        #             stop=2911
        #             id_length=len("SCMI01000001"),
        #     ]
        # ),
    ],
)
def test_can_resolve_wgs_ids_with_ena(accession, expected):
    assert wgs.resolve_ena_wgs(accession) == expected


# @pytest.mark.parametrize(
#     "given,accession,expected",
#     [
#         (
#             wgs.WgsSequence(
#                 wgs_id="LWDE02",
#                 kind=wgs.WgsSequenceKind.SEQUENCE,
#                 start=1,
#                 stop=23,
#                 id_length=len("LWDE02000001"),
#             ),
#             "LWDE01000000",
#             True,
#         ),
#         (
#             wgs.WgsSequence(
#                 wgs_id="LWDE02",
#                 kind=wgs.WgsSequenceKind.SEQUENCE,
#                 start=1,
#                 stop=23,
#                 id_length=len("LWDE02000001"),
#             ),
#             "LWDE01",
#             True,
#         ),
#         (
#             wgs.WgsSequence(
#                 wgs_id="ACTP02",
#                 kind=wgs.WgsSequenceKind.SEQUENCE,
#                 start=1,
#                 stop=23,
#                 id_length=len("ACTP02000001"),
#             ),
#             "LWDE01000000",
#             False,
#         ),
#     ],
# )
# def test_can_detect_if_within_one_version(given, accession, expected):
#     assert given.within_one_version(accession) == expected


@pytest.mark.parametrize(
    "accession,to_check,expected",
    [
        ("ALWZ04S0000000", "ALWZ04S0000000", True),
        ("ALWZ040000000", "ALWZ040000000", False),
        ("ALWZ04S0000001", "ALWZ04S", True),
        ("ALWZ04S", "ALWZ04S0000001", True),
        ("ALWZ04S", "ALWZ040000001", False),
        ("LKAM01", "LKAM01", True),
        ("LKAM01", "LKAM01000000", True),
        ("AMCR01", "LKAM01000000", False),
        ("WWJG01", "WWJG01002768.1", True),
        ("WWJG01", "WWJG01002768", True),
    ],
)
def test_can_correctly_match_wgs_ids(accession, to_check, expected):
    summary = wgs.resolve_ena_wgs(accession)
    assert summary.id_matches(to_check) == expected


@pytest.mark.parametrize(
    "accession,wgs_id,wgs_version,wgs_type",
    [
        ("ALWZ05", "ALWZ", 5, None),
        ("ACTP02", "ACTP", 2, None),
        ("ALWZ04S", "ALWZ", 4, "S"),
    ],
)
def test_wgs_summary_know_wgs_project_id(accession, wgs_id, wgs_version):
    summary = wgs.resolve_ena_wgs(accession)
    assert summary.wgs_accession == wgs_id
    assert summary.wgs_version == wgs_version
    assert summary.wgs_type == wgs_type
