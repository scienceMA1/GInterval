from ginterval import GInterval


def test_ginterval():
    i = GInterval(10, 20)
    assert len(i) == 10
    assert i.x == 10
    assert i.y == 20

    gi = GInterval(10, 50, chrom='chr1', name='GInterval1', strand='-')

    assert len(gi) == 40

    gi.set_strand_sensitivity(False)
    assert gi.start_position == 10
    assert gi.end_position == 49
    assert gi.start_index == 0
    assert gi.end_index == 39

    gi.set_strand_sensitivity(True)
    assert gi.start_position == 49
    assert gi.end_position == 10
    assert gi.start_index == 0
    assert gi.end_index == 39

    blocks = [(10, 20), (30, 35), (50, 65)]
    gi = GInterval(chrom='chr2', name='GInterval2', strand='-', blocks=blocks)

    assert gi.start_index == 0
    assert gi.end_index == 29
    assert len(gi) == 30

    assert gi.index(10) == 29
    gi.set_strand_sensitivity(False)
    assert gi.index(15) == 5
    gi.set_strand_sensitivity(True)
    assert gi.index(15) == 24

    assert gi.position(29) == 10

    intervals = [(10, 20), (30, 35), (50, 65)]
    thick = (15, 60)
    gi = GInterval(chrom='chr2', name='GInterval2', strand='-', blocks=intervals, thick=thick, rgb='read')

    assert gi.thick_start_position == 59
    assert gi.thick_end_position == 15
    assert gi.thick_start_index == 5
    assert gi.thick_end_index == 24
    assert gi.red is None

    g1 = gi[10:31]
    assert len(g1) == 11

    g1 = GInterval(blocks=[(10, 20), (20, 22), (100, 200)])
    g2 = GInterval(25, 27, chrom='chr2', name='GInterval2', strand='-', rgb='read')
    g3 = g1 + g2
    assert len(g3) == 114

    g1 = GInterval(blocks=[(10, 20), (30, 50), (90, 100)])
    g2 = GInterval(blocks=[(40, 50), (90, 95)])
    assert g1.coincide(g2)

    g1 = GInterval(blocks=[(10, 20), (30, 50), (90, 100)])
    g2 = GInterval(blocks=[(40, 52), (90, 95)])
    assert not g1.coincide(g2)
