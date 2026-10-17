from validation.release_certificate import all_required_pass

def test_certificate_requires_every_mandatory_pass():
    assert not all_required_pass({})
    assert all_required_pass({'one':{'status':'PASS'}})
    for status in ('FAIL','NOT_RUN','INCONCLUSIVE','OUT_OF_SCOPE'):
        assert not all_required_pass({'one':{'status':'PASS'},'two':{'status':status}})
