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
            start=1, stop=3033285)),
    (wgs.WgsSequenceKind.SEQUENCE, 'ALWZ040000001-ALWZ046221640',
        wgs.WgsInfo(wgs_id='ALWZ04', kind=wgs.WgsSequenceKind.SEQUENCE,
            start=1, stop=6221640)),
])
def test_wgs_info_parses_correctly(kind, raw, expected):
    info = wgs.WgsInfo.build(kind, raw)
    assert info == expected


@pytest.mark.parametrize('raw,expected', [
    ('ALWZ04S0000001', ('ALWZ04S', 1)),
    ('ALWZ04S6221640', ('ALWZ04S', 6221640)),
    ('ALWZ046221640', ('ALWZ04', 6221640)),
    ('LKAM01000000', ('LKAM01', 0)),
])
def test_can_get_endpoint_correctly(raw, expected):
    assert wgs.wgs_endpoint(raw) == expected


@pytest.mark.parametrize('info,expected', [
    (
        wgs.WgsInfo(wgs_id='ALWZ04', kind=wgs.WgsSequenceKind.SEQUENCE, start=1, stop=10),
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
        wgs.WgsInfo(wgs_id='ALWZ04S', kind=wgs.WgsSequenceKind.SCAFFOLD, start=1, stop=10),
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
])
def test_wgs_generates_expected_ids(info, expected):
    assert list(info.ids()) == expected


@pytest.mark.parametrize('accession,expected', [
    ('ALWZ05', [
        wgs.WgsInfo(wgs_id='ALWZ05', kind=wgs.WgsSequenceKind.SEQUENCE, start=1, stop=5154716),
        wgs.WgsInfo(wgs_id='ALWZ05S', kind=wgs.WgsSequenceKind.CONTIG, start=1, stop=2066983),
    ])
    
])
def test_can_resolve_wgs_ids_with_ena(accession, expected):
    assert wgs.resolve_ena_wgs(accession) == expected
