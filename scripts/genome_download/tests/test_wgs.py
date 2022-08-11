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

@pytest.mark.parametrize('kind,raw,expected', [
    (wgs.WgsSequenceKind.SCAFFOLD, 'ALWZ04S0000001-ALWZ04S3033285',
        wgs.WgsInfo(wgs_id='ALWZ04S', kind=wgs.WgsSequenceKind.SCAFFOLD,
            start=1, stop=3033285, id_length=len("ALWZ04S0000001"))),
    (wgs.WgsSequenceKind.SEQUENCE, 'ALWZ040000001-ALWZ046221640',
        wgs.WgsInfo(wgs_id='ALWZ04', kind=wgs.WgsSequenceKind.SEQUENCE,
            start=1, stop=6221640, id_length=len("ALWZ040000001"))),
    (wgs.WgsSequenceKind.SEQUENCE, 'AMCR01000001-AMCR01022363',
        wgs.WgsInfo(wgs_id='AMCR01', kind=wgs.WgsSequenceKind.SEQUENCE,
            start=1, stop=22363, id_length=len("AMCR01000001"))),
])
def test_wgs_info_parses_correctly(kind, raw, expected):
    info = wgs.WgsInfo.build(kind, raw)
    assert info == expected


@pytest.mark.parametrize('kind,raw', [
    (wgs.WgsSequenceKind.SCAFFOLD, 'ALWZ04S0000001-ALWZ04S3033285'),
    (wgs.WgsSequenceKind.SEQUENCE, 'ALWZ040000001-ALWZ046221640'),
    (wgs.WgsSequenceKind.SEQUENCE, 'AMCR01000001-AMCR01022363'),
])
def test_parsing_preseves_ids(kind, raw):
    info = wgs.WgsInfo.build(kind, raw)
    assert info.as_range() == raw


@pytest.mark.parametrize('info,id,expected', [
    (wgs.WgsInfo(wgs_id='ALWZ04S', kind=wgs.WgsSequenceKind.SCAFFOLD,
        start=1, stop=3033285, id_length=len('ALWZ04S0000001')),
        'ALWZ04S0000001', True,),
    (wgs.WgsInfo(wgs_id='ALWZ04S', kind=wgs.WgsSequenceKind.SCAFFOLD,
        start=1, stop=3033285, id_length=len('ALWZ04S0000001')),
        'ALWZ04S9999999', False,),
    (wgs.WgsInfo(wgs_id='AMCR01', kind=wgs.WgsSequenceKind.SEQUENCE, start=1, stop=22363,
        id_length=len("AMCR01000001")),
        "AMCR01000001", True,)
])
def test_wgs_knows_what_it_contains(info, id, expected):
    assert (id in info) == expected


@pytest.mark.parametrize('raw,expected', [
    ('ALWZ04S0000001', ('ALWZ04S', 1)),
    ('ALWZ04S6221640', ('ALWZ04S', 6221640)),
    ('ALWZ046221640', ('ALWZ04', 6221640)),
    ('LKAM01000000', ('LKAM01', 0)),
    ('AMCR01000000', ('AMCR01', 0)),
])
def test_can_get_endpoint_correctly(raw, expected):
    assert wgs.wgs_endpoint(raw) == expected


@pytest.mark.parametrize('info,expected', [
    (
        wgs.WgsInfo(wgs_id='ALWZ04', kind=wgs.WgsSequenceKind.SEQUENCE, start=1, stop=10,
            id_length=len('ALWZ040000001')),
        [
            'ALWZ040000001',
            'ALWZ040000002',
            'ALWZ040000003',
            'ALWZ040000004',
            'ALWZ040000005',
            'ALWZ040000006',
            'ALWZ040000007',
            'ALWZ040000008',
            'ALWZ040000009',
            'ALWZ040000010',
        ]
    ),
    (
        wgs.WgsInfo(wgs_id='ALWZ04S', kind=wgs.WgsSequenceKind.SCAFFOLD, start=1, stop=10,
            id_length=len('ALWZ04S0000001')),
        [
            'ALWZ04S0000001',
            'ALWZ04S0000002',
            'ALWZ04S0000003',
            'ALWZ04S0000004',
            'ALWZ04S0000005',
            'ALWZ04S0000006',
            'ALWZ04S0000007',
            'ALWZ04S0000008',
            'ALWZ04S0000009',
            'ALWZ04S0000010',
        ]
    ),
    (
        wgs.WgsInfo(wgs_id='AMCR01', kind=wgs.WgsSequenceKind.SCAFFOLD, start=1, stop=10,
            id_length=len('AMCR010000001')),
        [
            'AMCR010000001',
            'AMCR010000002',
            'AMCR010000003',
            'AMCR010000004',
            'AMCR010000005',
            'AMCR010000006',
            'AMCR010000007',
            'AMCR010000008',
            'AMCR010000009',
            'AMCR010000010',
        ]
    ),
])
def test_wgs_generates_expected_ids(info, expected):
    assert list(info.ids()) == expected


@pytest.mark.parametrize('accession,expected', [
    ('ALWZ05', [
        wgs.WgsInfo(wgs_id='ALWZ05', kind=wgs.WgsSequenceKind.SEQUENCE, start=1, stop=5154716,
            id_length=13),
    ])
    
])
def test_can_resolve_wgs_ids_with_ena(accession, expected):
    assert wgs.resolve_ena_wgs(accession) == expected
