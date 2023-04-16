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
    "raw,expected",
    [
        (
            "JABDTM02",
            wgs.WgsPrefix(wgs_id="JABDTM", wgs_version=2, scaffold=False, length=8),
        ),
        (
            "ACTP02",
            wgs.WgsPrefix(wgs_id="ACTP", wgs_version=2, scaffold=False, length=6),
        ),
        (
            "ALWZ04S",
            wgs.WgsPrefix(wgs_id="ALWZ", wgs_version=4, scaffold=True, length=7),
        ),
        (
            "AAFI02000000",
            wgs.WgsPrefix(wgs_id="AAFI", wgs_version=2, scaffold=False, length=6),
        ),
    ],
)
def test_can_build_wgs_prefix_correctly(raw, expected):
    assert wgs.WgsPrefix.build(raw) == expected


@pytest.mark.parametrize(
    "raw",
    [
        "JABDTM02",
        "ACTP02",
        "ALWZ04S",
    ],
)
def test_wgs_prefix_can_roundtrip(raw):
    assert wgs.WgsPrefix.build(raw).to_wgs_string() == raw


@pytest.mark.parametrize(
    "raw,expected",
    [
        (
            "AAFI02000000",
            "AAFI02",
        ),
    ],
)
def test_truncates_long_accession_in_roundtrip(raw, expected):
    assert wgs.WgsPrefix.build(raw).to_wgs_string() == expected


@pytest.mark.parametrize(
    "first,second,expected",
    [
        ("JABDTM02", "JABDTM02", True),
        ("JABDTM02", "ACTP02", False),
        ("ACTP02", "JABDTM02", False),
    ],
)
def test_wgs_prefix_can_detect_if_they_match(first, second, expected):
    left = wgs.WgsPrefix.build(first)
    right = wgs.WgsPrefix.build(second)
    assert left.matches(right) == expected


@pytest.mark.parametrize(
    "first,second,expected",
    [
        ("JABDTM02", "JABDTM02", True),
        ("JABDTM02", "ACTP02", False),
        ("ACTP02", "JABDTM02", False),
        ("JABDTM01", "JABDTM02", True),
        ("JABDTM02", "JABDTM01", True),
    ],
)
def test_wgs_prefix_can_detect_if_they_match_within_a_version(first, second, expected):
    left = wgs.WgsPrefix.build(first)
    right = wgs.WgsPrefix.build(second)
    assert left.matches(right, within_one_version=True) == expected


@pytest.mark.parametrize(
    "raw,expected",
    [
        (
            "ALWZ04S3033285",
            wgs.WgsSequenceId(
                prefix=wgs.WgsPrefix(
                    wgs_id="ALWZ", wgs_version=4, scaffold=True, length=7
                ),
                sequence_index=3033285,
                sequence_version=None,
                length=14,
            ),
        ),
        (
            "BAMX01000058.1",
            wgs.WgsSequenceId(
                prefix=wgs.WgsPrefix(
                    wgs_id="BAMX",
                    wgs_version=1,
                    scaffold=False,
                    length=6,
                ),
                sequence_index=58,
                sequence_version="1",
                length=14,
            ),
        ),
    ],
)
def test_can_parse_wgs_sequence_id(raw, expected):
    assert wgs.WgsSequenceId.build(raw) == expected


@pytest.mark.parametrize(
    "raw",
    [
        ("ALWZ04S3033285"),
        ("BAMX01000058.1"),
    ],
)
def test_can_roundtrip_sequence_ids(raw):
    assert wgs.WgsSequenceId.build(raw).to_wgs_string() == raw


@pytest.mark.parametrize(
    "raw,expected",
    [
        (
            "ALWZ04S0000001-ALWZ04S3033285",
            wgs.WgsSequenceRange(
                prefix=wgs.WgsPrefix.build("ALWZ04S"),
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=wgs.WgsSequenceId.build("ALWZ04S0000001"),
                stop=wgs.WgsSequenceId.build("ALWZ04S3033285"),
            ),
        ),
        (
            "ALWZ040000001-ALWZ046221640",
            wgs.WgsSequenceRange(
                prefix=wgs.WgsPrefix.build("ALWZ04"),
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=wgs.WgsSequenceId.build("ALWZ040000001"),
                stop=wgs.WgsSequenceId.build("ALWZ046221640"),
            ),
        ),
        (
            "AMCR01000001-AMCR01022363",
            wgs.WgsSequenceRange(
                prefix=wgs.WgsPrefix.build("AMCR01"),
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=wgs.WgsSequenceId.build("AMCR01000001"),
                stop=wgs.WgsSequenceId.build("AMCR01022363"),
            ),
        ),
        (
            "ACTP02000001-ACTP02000023",
            wgs.WgsSequenceRange(
                prefix=wgs.WgsPrefix.build("ACTP02"),
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=wgs.WgsSequenceId.build("ACTP02000001"),
                stop=wgs.WgsSequenceId.build("ACTP02000023"),
            ),
        ),
        (
            "JABDTM020000001-JABDTM020000010",
            wgs.WgsSequenceRange(
                prefix=wgs.WgsPrefix.build("JABDTM02"),
                kind=wgs.WgsSequenceKind.SEQUENCE,
                start=wgs.WgsSequenceId.build("JABDTM020000001"),
                stop=wgs.WgsSequenceId.build("JABDTM020000010"),
            ),
        ),
    ],
)
def test_can_parse_wgs_sequence_set(raw, expected):
    assert wgs.WgsSequenceRange.build(wgs.WgsSequenceKind.SEQUENCE, raw) == expected


@pytest.mark.parametrize(
    "raw",
    [
        ("ALWZ04S0000001-ALWZ04S3033285"),
        ("ALWZ040000001-ALWZ046221640"),
        ("AMCR01000001-AMCR01022363"),
        ("ACTP02000001-ACTP02000023"),
        ("JABDTM020000001-JABDTM020000010"),
    ],
)
def test_can_roundtrip_wgs_sequence_sets(raw):
    assert (
        wgs.WgsSequenceRange.build(wgs.WgsSequenceKind.SEQUENCE, raw).to_wgs_string()
        == raw
    )


@pytest.mark.parametrize(
    "raw,id,expected",
    [
        (
            "ALWZ04S0000001-ALWZ04S3033285",
            "ALWZ04S0000001",
            True,
        ),
        (
            "ALWZ04S0000001-ALWZ04S3033285",
            "ALWZ040000001",
            False,
        ),
        (
            "AMCR01000001-AMCR01022363",
            "JABDTM020000011",
            False,
        ),
        ("ACTP02000001-ACTP02000023", "ACTP02000020", True),
        (
            "JABDTM020000001-JABDTM020000010",
            "JABDTM020000001",
            True,
        ),
        (
            "JABDTM020000001-JABDTM020000010",
            "JABDTM020000010",
            True,
        ),
        (
            "JABDTM020000001-JABDTM020000010",
            "JABDTM020000011",
            False,
        ),
    ],
)
def test_can_tell_if_sequence_id_is_member(raw, id, expected):
    assert (
        wgs.WgsSequenceRange.build(wgs.WgsSequenceKind.SEQUENCE, raw).includes(id)
        == expected
    )


@pytest.mark.parametrize(
    "raw,expected",
    [
        (
            "ALWZ04S0000001",
            True,
        ),
        (
            "ALWZ04S0000001-ALWZ04S0000001",
            True,
        ),
        (
            "ALWZ04S0000001-ALWZ04S3033285",
            False,
        ),
    ],
)
def test_can_tell_if_sequence_range_is_single(raw, expected):
    range = wgs.WgsSequenceRange.build(wgs.WgsSequenceKind.SEQUENCE, raw)
    assert range.is_single_range() == expected


# @pytest.mark.parametrize(
#     "kind,raw,expected",
#     [
#         (
#             wgs.WgsSequenceKind.SCAFFOLD,
#             "ALWZ04S0000001-ALWZ04S3033285",
#             "ALWZ",
#             "04S",
#         ),
#         (
#             wgs.WgsSequenceKind.SEQUENCE,
#             "ALWZ040000001-ALWZ046221640",
#             "ALWZ",
#             "04",
#         ),
#         (
#             wgs.WgsSequenceKind.SEQUENCE,
#             "AMCR01000001-AMCR01022363",
#             "AMCR",
#             "01",
#         ),
#         (
#             wgs.WgsSequenceKind.CONTIG,
#             "ACTP02000001-ACTP02000023",
#             "ACTP",
#             "02",
#         ),
#         (
#             wgs.WgsSequenceKind.SEQUENCE,
#             "JABDTM020000001-JABDTM020000010",
#             "JABDTM",
#             "02",
#         )
#     ],
# )
# def test_can_correctly_fetch_wgs_info(kind, raw, prefix, version):
#     info = wgs.WgsSequenceSet.build(kind, raw)
#     assert info.wgs_prefix == prefix
#     assert info.wgs_version == version
#
#
#
# # @pytest.mark.parametrize(
# #     "raw,expected",
# #     [
# #         ("KE150405-KE150411", wgs.ContigInfo(prefix="KE", start=150405, stop=150411)),
# #     ],
# # )
# # def test_can_parse_scaffolds(raw, expected):
# #     info = wgs.ContigInfo.build(raw)
# #     assert info == expected
# #
# #
#
#
# # @pytest.mark.parametrize(
# #     "raw,expected",
# #     [
# #         ("KE150405-KE150411", wgs.ContigInfo(prefix="KE", start=150405, stop=150411)),
# #     ],
# # )
# # def test_can_parse_scaffolds(raw, expected):
# #     info = wgs.ContigInfo.build(raw)
# #     assert info == expected
#
#
# # @pytest.mark.parametrize(
# #     "kind,raw",
# #     [
# #         (wgs.WgsSequenceKind.SCAFFOLD, "ALWZ04S0000001-ALWZ04S3033285"),
# #         (wgs.WgsSequenceKind.SEQUENCE, "ALWZ040000001-ALWZ046221640"),
# #         (wgs.WgsSequenceKind.SEQUENCE, "AMCR01000001-AMCR01022363"),
# #     ],
# # )
# # def test_parsing_preseves_ids(kind, raw):
# #     info = wgs.WgsSequence.build(kind, raw)
# #     assert info.as_range() == raw
#
#
# # @pytest.mark.parametrize(
# #     "info,id,expected",
# #     [
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="ALWZ04S",
# #                 kind=wgs.WgsSequenceKind.SCAFFOLD,
# #                 start=1,
# #                 stop=3033285,
# #                 id_length=len("ALWZ04S0000001"),
# #             ),
# #             "ALWZ04S0000001",
# #             True,
# #         ),
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="ALWZ04S",
# #                 kind=wgs.WgsSequenceKind.SCAFFOLD,
# #                 start=1,
# #                 stop=3033285,
# #                 id_length=len("ALWZ04S0000001"),
# #             ),
# #             "ALWZ04S9999999",
# #             False,
# #         ),
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="AMCR01",
# #                 kind=wgs.WgsSequenceKind.SEQUENCE,
# #                 start=1,
# #                 stop=22363,
# #                 id_length=len("AMCR01000001"),
# #             ),
# #             "AMCR01000001",
# #             True,
# #         ),
# #     ],
# # )
# # def test_wgs_knows_what_it_contains(info, id, expected):
# #     assert (id in info) == expected
#
#
# @pytest.mark.parametrize(
#     "raw,expected",
#     [
#         ("ALWZ04S0000001", wgs.WgsSequenceId(
#             wgs_id="ALWZ",
#             wgs_version="04S",
#             sequence_index=1,
#             sequence_version=None,
#         )),
#         ("ALWZ04S6221640", wgs.WgsSequenceId(
#             wgs_id="ALWZ",
#             wgs_version="04S",
#             sequence_index=6221640,
#             sequence_version=None,
#         )),
#         ("ALWZ046221640", wgs.WgsSequenceId(
#             wgs_id="ALWZ",
#             wgs_version="04",
#             sequence_index=6221640,
#             sequence_version=None,
#         )),
#         ("LKAM01000000", wgs.WgsSequenceId(
#             wgs_id="LKAM",
#             wgs_version="01",
#             sequence_index=0,
#             sequence_version=None,
#         )),
#         ("AMCR01000000", wgs.WgsSequenceId(
#             wgs_id="AMCR",
#             wgs_version="01",
#             sequence_index=0,
#             sequence_version=None,
#         )),
#         ("BAMX01000058.1", wgs.WgsSequenceId(
#             wgs_id="BAMX",
#             wgs_version="01",
#             sequence_index=58,
#             sequence_version='1',
#         )),
#     ],
# )
# def test_can_get_endpoint_correctly(raw, expected):
#     assert wgs.wgs_endpoint(raw) == expected
#
#
# # @pytest.mark.parametrize(
# #     "info,expected",
# #     [
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="ALWZ04",
# #                 kind=wgs.WgsSequenceKind.SEQUENCE,
# #                 start=1,
# #                 stop=10,
# #                 id_length=len("ALWZ040000001"),
# #             ),
# #             [
# #                 "ALWZ040000001",
# #                 "ALWZ040000002",
# #                 "ALWZ040000003",
# #                 "ALWZ040000004",
# #                 "ALWZ040000005",
# #                 "ALWZ040000006",
# #                 "ALWZ040000007",
# #                 "ALWZ040000008",
# #                 "ALWZ040000009",
# #                 "ALWZ040000010",
# #             ],
# #         ),
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="ALWZ04S",
# #                 kind=wgs.WgsSequenceKind.SCAFFOLD,
# #                 start=1,
# #                 stop=10,
# #                 id_length=len("ALWZ04S0000001"),
# #             ),
# #             [
# #                 "ALWZ04S0000001",
# #                 "ALWZ04S0000002",
# #                 "ALWZ04S0000003",
# #                 "ALWZ04S0000004",
# #                 "ALWZ04S0000005",
# #                 "ALWZ04S0000006",
# #                 "ALWZ04S0000007",
# #                 "ALWZ04S0000008",
# #                 "ALWZ04S0000009",
# #                 "ALWZ04S0000010",
# #             ],
# #         ),
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="AMCR01",
# #                 kind=wgs.WgsSequenceKind.SCAFFOLD,
# #                 start=1,
# #                 stop=10,
# #                 id_length=len("AMCR010000001"),
# #             ),
# #             [
# #                 "AMCR010000001",
# #                 "AMCR010000002",
# #                 "AMCR010000003",
# #                 "AMCR010000004",
# #                 "AMCR010000005",
# #                 "AMCR010000006",
# #                 "AMCR010000007",
# #                 "AMCR010000008",
# #                 "AMCR010000009",
# #                 "AMCR010000010",
# #             ],
# #         ),
# #     ],
# # )
# # def test_wgs_generates_expected_ids(info, expected):
# #     assert list(info.ids()) == expected
#
#
# # @pytest.mark.parametrize(
# #     "accession,expected",
# #     [
# #         (
# #             "ALWZ05",
# #             (
# #                 [],
# #                 [
# #                     wgs.WgsSequence(
# #                         wgs_id="ALWZ05",
# #                         kind=wgs.WgsSequenceKind.SEQUENCE,
# #                         start=1,
# #                         stop=5154716,
# #                         id_length=13,
# #                     ),
# #                     wgs.WgsSequence(
# #                         wgs_id="ALWZ05S",
# #                         kind=wgs.WgsSequenceKind.CONTIG,
# #                         start=1,
# #                         stop=2066983,
# #                         id_length=len("ALWZ05S0000001"),
# #                     ),
# #                 ],
# #             ),
# #         ),
# #         (
# #             "ACTP02",
# #             (
# #                 [wgs.ContigInfo(prefix="KE", start=150405, stop=150411)],
# #                 [
# #                     wgs.WgsSequence(
# #                         wgs_id="ACTP02",
# #                         kind=wgs.WgsSequenceKind.SEQUENCE,
# #                         start=1,
# #                         stop=23,
# #                         id_length=len("ACTP02000001"),
# #                     ),
# #                 ],
# #             ),
# #         ),
# #         # (
# #         #     "SCMI01",
# #         #     [
# #         #         wgs.WgsSequence(
# #         #             wgs_id='SCMI01',
# #         #             kind=wgs.WgsSequenceKind.SEQUENCE,
# #         #             start=1,
# #         #             stop=2911
# #         #             id_length=len("SCMI01000001"),
# #         #     ]
# #         # ),
# #     ],
# # )
# # def test_can_resolve_wgs_ids_with_ena(accession, expected):
# #     assert wgs.resolve_ena_wgs(accession) == expected
#
#
# @pytest.mark.parametrize(
#     "accession,to_check,expected",
#     [
#         ("ALWZ04S0000000", "ALWZ04S0000000", True),
#         ("ALWZ040000000", "ALWZ040000000", False),
#         ("ALWZ04S0000001", "ALWZ04S", True),
#         ("ALWZ04S", "ALWZ04S0000001", True),
#         ("ALWZ04S", "ALWZ040000001", False),
#         ("LKAM01", "LKAM01000000", True),
#         ("AMCR01", "LKAM01000000", False),
#         ("WWJG01", "WWJG01002768.1", True),
#         ("WWJG01", "WWJG01002768", True),
#     ],
# )
# def test_can_correctly_match_wgs_ids(accession, to_check, expected):
#     (contigs, ena_info) = wgs.resolve_ena_wgs(accession)
#     wgs_id = ena_info[0].wgs_id
#     endpoint = wgs.wgs_endpoint(to_check)
#     summary = wgs.WgsSummary(
#         wgs_id=wgs_id,
#         contigs=contigs,
#         sequences=ena_info,
#         ncbi_ids=[],
#     )
#     assert summary.id_matches(endpoint) == expected
#
#
# # @pytest.mark.skip()
# # @pytest.mark.parametrize(
# #     "accession,wgs_id,wgs_version,wgs_type",
# #     [
# #         ("ALWZ05", "ALWZ", 5, None),
# #         ("ACTP02", "ACTP", 2, None),
# #         ("ALWZ04S", "ALWZ", 4, "S"),
# #     ],
# # )
# # def test_wgs_summary_know_wgs_project_id(accession, wgs_id, wgs_version):
# #     summary = wgs.resolve_ena_wgs(accession)
# #     assert summary.wgs_accession == wgs_id
# #     assert summary.wgs_version == wgs_version
# #     assert summary.wgs_type == wgs_type
#
#
# # @pytest.mark.parametrize(
# #     "given,accession,expected",
# #     [
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="LWDE02",
# #                 kind=wgs.WgsSequenceKind.SEQUENCE,
# #                 start=1,
# #                 stop=23,
# #                 id_length=len("LWDE02000001"),
# #             ),
# #             "LWDE01000000",
# #             True,
# #         ),
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="LWDE02",
# #                 kind=wgs.WgsSequenceKind.SEQUENCE,
# #                 start=1,
# #                 stop=23,
# #                 id_length=len("LWDE02000001"),
# #             ),
# #             "LWDE01",
# #             True,
# #         ),
# #         (
# #             wgs.WgsSequence(
# #                 wgs_id="ACTP02",
# #                 kind=wgs.WgsSequenceKind.SEQUENCE,
# #                 start=1,
# #                 stop=23,
# #                 id_length=len("ACTP02000001"),
# #             ),
# #             "LWDE01000000",
# #             False,
# #         ),
# #     ],
# # )
# # def test_can_detect_if_within_one_version(given, accession, expected):
# #     assert given.within_one_version(accession) == expected
