# -*- coding: utf-8 -*-
"""Headless GUI smoke test for single-stranded loopout rendering in PathView.

Run with:
    pytest -c cadnano/tests/pytestgui.ini cadnano/tests/loopouttest_gui.py

Model imports are deliberately deferred until after the Qt application boots
so the PyQt signal proxy is configured before CNObject subclasses are defined.
A single test boots one QApplication (Qt only supports one per process).
"""
import pytest
from PyQt5.QtCore import Qt

from cnguitestcase import GUITestApp
from cadnano.proxies.cnenum import ViewSendEnum


@pytest.fixture()
def cnapp():
    app = GUITestApp()
    yield app
    app.tearDown()


def _path_part_item(cnapp, part):
    root = cnapp.app_window.views[ViewSendEnum.PATH].rootItem()
    part_items = [item for item in root.partItems()
                  if item.part() is part]
    assert part_items
    return part_items[0]


def _strand_item(vhi, model_strand):
    for item in vhi.getStrandItems():
        if item.strand() is model_strand:
            return item
    return None


def test_loopout_rendered_and_updated_in_pathview(cnapp):
    # import after GUI bootstrap so proxyConfigure('PyQt') has run
    from nucleicacidparttest import create3Helix
    part = create3Helix(cnapp.document, [0, 0, 1], 42)
    fwd0, _ = part.getStrandSets(0)
    fwd1, _ = part.getStrandSets(1)

    # 1. plain crossover: no loopout badge, solid pen
    s5p = fwd0.createStrand(0, 10)
    s3p = fwd1.createStrand(20, 30)
    part.createXover(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime())
    cnapp.processEvents()

    part_item = _path_part_item(cnapp, part)
    vhi0 = part_item.vhItemForIdNum(0)
    si = _strand_item(vhi0, s5p)
    assert si is not None
    xover_item = si.xover_3p_end
    assert xover_item._strand5p is s5p
    assert xover_item._loopout_label is None or \
        not xover_item._loopout_label.isVisible()
    assert xover_item.pen().style() == Qt.SolidLine

    # 2. crossover -> loopout: dashed pen and editable length badge appear
    part.setLoopoutLength(s5p, 3)
    cnapp.processEvents()
    assert xover_item._loopout_label is not None
    assert xover_item._loopout_label.isVisible()
    assert xover_item._loopout_label.text() == "L3"
    assert xover_item.pen().style() == Qt.DashLine

    # 3. length edits flow from the model back into the view
    part.setLoopoutLength(s5p, 5)
    cnapp.processEvents()
    assert xover_item._loopout_label.text() == "L5"
    assert s5p.loopoutLength() == 5

    # 4. loopout -> crossover: solid pen returns, badge hidden
    part.setLoopoutLength(s5p, 0)
    cnapp.processEvents()
    assert not s5p.hasLoopout()
    assert xover_item.pen().style() == Qt.SolidLine
    assert not xover_item._loopout_label.isVisible()

    # 5. a freshly created loopout renders the same way
    part.setLoopoutLength(s5p, 2)
    cnapp.processEvents()
    assert xover_item._loopout_label.text() == "L2"
    assert xover_item.pen().style() == Qt.DashLine
