# -*- coding: utf-8 -*-
"""Tests for single-stranded loopout support.

A loopout is a connection between the free 3' end of a bound strand and the
free 5' end of another bound strand carrying one or more unpaired bases that
occupy no virtual helix grid point. Length 0 is a plain crossover.
"""
import io
import json

import pytest

from cntestcase import cnapp
from nucleicacidparttest import create3Helix

from cadnano.fileio import encode as cnencode
from cadnano.fileio import v3decode
from cadnano.fileio.v2encode import LoopoutLegacyExportError


def _two_fwd_strands(part, length=42):
    """Two forward strands on helices 0 and 1 with free 3' and 5' ends."""
    fwd0, _ = part.getStrandSets(0)
    fwd1, _ = part.getStrandSets(1)
    s5p = fwd0.createStrand(0, 10)
    s3p = fwd1.createStrand(20, 30)
    return s5p, s3p


def test_loopout_creation_and_validation(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    us = part.undoStack()

    # illegal lengths are rejected and leave no undo record
    idx_before = len(us.undostack)
    assert part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), -1) is False
    assert part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 1.5) is False
    assert part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), True) is False
    assert len(us.undostack) == idx_before
    assert s5p.connection3p() is None and s3p.connection5p() is None

    # endpoints that are not the actual free bound ends are rejected
    assert part.createLoopout(s5p, s5p.lowIdx(), s3p, s3p.idx5Prime(), 3) is False
    assert part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.highIdx(), 3) is False
    assert len(us.undostack) == idx_before

    # circularization of an oligo is rejected
    assert part.createLoopout(s5p, s5p.idx3Prime(), s5p, s5p.idx5Prime(), 3) is False
    assert len(us.undostack) == idx_before

    # a valid loopout is created and counts in the oligo length
    assert part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 3) is True
    assert s5p.connection3p() is s3p
    assert s3p.connection5p() is s5p
    assert s5p.hasLoopout() and s5p.loopoutLength() == 3
    # bound geometry is untouched
    assert s5p.idxs() == (0, 10) and s3p.idxs() == (20, 30)
    merged = s5p.oligo()
    assert merged is s3p.oligo()
    assert merged.length() == s5p.totalLength() + 3 + s3p.totalLength()

    # an endpoint now occupied by another connection is rejected
    other5p = part.getStrandSets(2)[0].createStrand(0, 5)
    assert part.createLoopout(other5p, other5p.idx3Prime(),
                               s3p, s3p.idx5Prime(), 2) is False
    assert other5p.connection3p() is None


def test_loopout_creation_preserves_side_sequences(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    # assign sequences to the two separate oligos first
    seq5 = 'ACGTACGTACA'[:s5p.totalLength()]
    seq3 = 'GTACGTACGTG'[:s3p.totalLength()]
    s5p.oligo().applySequence(seq5)
    s3p.oligo().applySequence(seq3)
    assert s5p.sequence() == seq5 and s3p.sequence() == seq3

    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 2)
    # both bound fragments keep their assigned bases; the loopout is unknown
    assert s5p.sequence() == seq5 and s3p.sequence() == seq3
    assert s5p.loopoutSequence() == '  '
    merged_seq = s5p.oligo().sequence()
    assert merged_seq.startswith(seq5) and seq3 in merged_seq


def test_loopout_undo_redo(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    oligo3p = s3p.oligo()
    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 4)

    part.undoStack().undo()
    assert s5p.connection3p() is None and s3p.connection5p() is None
    assert not s5p.hasLoopout()
    assert s5p.oligo() is not s3p.oligo()
    assert s3p.oligo() is oligo3p

    part.undoStack().redo()
    assert s5p.connection3p() is s3p and s5p.loopoutLength() == 4
    assert s5p.oligo() is s3p.oligo()
    assert s5p.oligo().length() == s5p.totalLength() + 4 + s3p.totalLength()


def test_loopout_length_change_preserves_bases(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 5)
    oligo = s5p.oligo()
    full_len = oligo.length()

    # assign a sequence long enough to cover bound bases plus the loopout
    bound5 = s5p.totalLength()
    seq = 'ACGT' * ((full_len // 4) + 1)
    oligo.applySequence(seq[:full_len])
    expected_loop = seq[bound5:bound5 + 5]
    assert s5p.loopoutSequence() == expected_loop
    bound_seq_5 = s5p.sequence()
    bound_seq_3 = s3p.sequence()

    # shrink: 5' prefix is retained, side sequences untouched
    assert part.setLoopoutLength(s5p, 2) is True
    assert s5p.loopoutSequence() == expected_loop[:2]
    assert s5p.sequence() == bound_seq_5 and s3p.sequence() == bound_seq_3
    assert oligo.length() == full_len - 3

    # grow: new positions are unknown
    assert part.setLoopoutLength(s5p, 4) is True
    assert s5p.loopoutSequence() == expected_loop[:2] + '  '
    assert oligo.length() == full_len - 1

    # invalid changes rejected without an undo record
    us = part.undoStack()
    idx = len(us.undostack)
    assert part.setLoopoutLength(s5p, -2) is False
    assert part.setLoopoutLength(s5p, 1.0) is False
    assert len(us.undostack) == idx

    # one undo per operation restores the complete prior state
    us.undo()
    assert s5p.loopoutLength() == 2 and s5p.loopoutSequence() == expected_loop[:2]
    us.undo()
    assert s5p.loopoutLength() == 5 and s5p.loopoutSequence() == expected_loop


def test_crossover_loopout_interconversion(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    len5, len3 = s5p.totalLength(), s3p.totalLength()

    # length 0 creates a plain crossover through the xover machinery
    assert part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 0) is True
    assert s5p.connection3p() is s3p and not s5p.hasLoopout()
    assert s5p.oligo().length() == len5 + len3

    # crossover -> loopout
    assert part.setLoopoutLength(s5p, 3) is True
    assert s5p.hasLoopout() and s5p.loopoutLength() == 3
    assert s5p.oligo().length() == len5 + len3 + 3

    # loopout -> crossover (length 0), connection retained
    assert part.setLoopoutLength(s5p, 0) is True
    assert s5p.connection3p() is s3p and not s5p.hasLoopout()
    assert s5p.oligo().length() == len5 + len3

    part.undoStack().undo()
    assert s5p.loopoutLength() == 3
    part.undoStack().undo()
    assert not s5p.hasLoopout() and s5p.connection3p() is s3p


def test_oligo_sequence_and_csv_export(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 4)
    oligo = s5p.oligo()
    total = oligo.length()
    seq = ('ACGT' * ((total // 4) + 1))[:total]
    oligo.applySequence(seq)

    # 5'->3' oligo sequence includes the loopout in the right position
    assert oligo.sequence() == seq
    b5 = s5p.totalLength()
    assert s5p.loopoutSequence() == seq[b5:b5 + 4]

    # CSV export marks unknown positions as '?'; here everything is known
    csv_text = part.getSequences()
    assert seq in csv_text

    # a loopout with unknown bases exports '?'
    part.setLoopoutLength(s5p, 6)  # two freshly added unknown positions
    csv_text = part.getSequences()
    assert '?' in csv_text
    # the unknown loopout positions sit between the two bound fragments
    row = [line for line in csv_text.splitlines() if seq[:b5] in line][0]
    exported_seq = row.split(',')[4]
    assert exported_seq[b5:b5 + 4] == seq[b5:b5 + 4]
    assert exported_seq[b5 + 4:b5 + 6] == '??'


def test_auto_complement_skips_loopout(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    # overlapping duplex on helix 0, plus a loopout leaving its 3' end
    fwd0, rev0 = part.getStrandSets(0)
    fwd1, _ = part.getStrandSets(1)
    s5p = fwd0.createStrand(0, 10)
    comp = rev0.createStrand(0, 10)
    s3p = fwd1.createStrand(20, 25)
    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 3)

    oligo = s5p.oligo()
    total = oligo.length()
    oligo.applySequence(('ACGT' * (total // 4 + 1))[:total])

    # the complement strand only receives paired bound bases (11), no loopout
    assert comp.oligo().length() == comp.totalLength()
    assert len(comp.sequence()) == comp.totalLength()


def test_multiple_loopouts_reverse_and_insertions(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    fwd0, rev0 = part.getStrandSets(0)
    fwd1, rev1 = part.getStrandSets(1)
    fwd2, _ = part.getStrandSets(2)

    a = fwd0.createStrand(0, 5)
    b = rev1.createStrand(10, 15)      # reverse fragment
    c = fwd2.createStrand(0, 4)
    part.createLoopout(a, a.idx3Prime(), b, b.idx5Prime(), 2)
    part.createLoopout(b, b.idx3Prime(), c, c.idx5Prime(), 5)
    # mix in an insertion and a skip
    a.addInsertion(2, 3)
    b.addInsertion(12, -1)

    oligo = a.oligo()
    assert b.oligo() is oligo and c.oligo() is oligo
    expected = (a.totalLength() + 2 + b.totalLength() + 5 + c.totalLength())
    assert oligo.length() == expected

    total = oligo.length()
    seq = ('ACGT' * (total // 4 + 1))[:total]
    oligo.applySequence(seq)
    assert oligo.sequence() == seq
    # loopout bases land at the right offsets
    off_a = a.totalLength()
    assert a.loopoutSequence() == seq[off_a:off_a + 2]
    off_b = off_a + 2 + b.totalLength()
    assert b.loopoutSequence() == seq[off_b:off_b + 5]

    # export works and every character position is accounted for
    csv_text = part.getSequences()
    exported = [line for line in csv_text.splitlines()
                if line.startswith('%d[' % a.idNum())][0].split(',')[4]
    assert len(exported) == total


def test_split_keeps_loopout_on_3p_fragment(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 3)
    oligo = s5p.oligo()
    length_before = oligo.length()

    # split the 5' strand well inside its range
    fwd0 = part.getStrandSets(0)[0]
    assert fwd0.splitStrand(s5p, 4) is True
    low = fwd0.getStrand(0)
    high = fwd0.getStrand(10)
    # the 3' fragment inherits the loopout, the 5' fragment is free
    assert high.connection3p() is s3p and high.loopoutLength() == 3
    assert low.connection3p() is None and not low.hasLoopout()
    # no connection was lost and total bases are conserved
    assert s3p.connection5p() is high
    assert high.oligo().length() + low.oligo().length() == length_before

    part.undoStack().undo()
    assert fwd0.getStrand(0) is s5p
    assert s5p.connection3p() is s3p and s5p.loopoutLength() == 3
    assert s5p.oligo().length() == length_before


def test_merge_and_remove_strand_cleanup(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    fwd0, _ = part.getStrandSets(0)
    fwd1, _ = part.getStrandSets(1)
    a = fwd0.createStrand(0, 5)
    mid = fwd0.createStrand(6, 10)
    s3p = fwd1.createStrand(20, 25)
    part.createLoopout(mid, mid.idx3Prime(), s3p, s3p.idx5Prime(), 4)
    oligo_len = mid.oligo().length()
    a_len = a.oligo().length()

    # merge the adjacent fragments; the loopout follows the 3' fragment
    assert fwd0.mergeStrands(a, mid) is True
    merged = fwd0.getStrand(0)
    assert merged.connection3p() is s3p and merged.loopoutLength() == 4
    # merging joins the two oligos, so the length grows by a's contribution
    assert merged.oligo().length() == oligo_len + a_len
    part.undoStack().undo()
    assert mid.loopoutLength() == 4 and mid.connection3p() is s3p
    assert mid.oligo().length() == oligo_len

    # removing the 3' strand of a loopout leaves no dangling connection
    part.undoStack().redo()
    merged = fwd0.getStrand(0)
    merged_len = merged.oligo().length()
    fwd1.removeStrand(s3p)
    assert merged.connection3p() is None and not merged.hasLoopout()
    part.undoStack().undo()
    assert merged.connection3p() is s3p and merged.loopoutLength() == 4
    assert s3p.connection5p() is merged
    assert merged.oligo().length() == merged_len


def test_remove_connection_restores_state(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 3)
    oligo5p = s5p.oligo()
    len5p = oligo5p.length()

    part.removeXover(s5p, s3p)  # removes any connection, loopout included
    assert s5p.connection3p() is None and s3p.connection5p() is None
    assert not s5p.hasLoopout()
    assert s5p.oligo() is not s3p.oligo()
    assert s5p.oligo().length() + s3p.oligo().length() == len5p

    part.undoStack().undo()
    assert s5p.connection3p() is s3p and s5p.loopoutLength() == 3
    assert s5p.oligo() is oligo5p and s5p.oligo().length() == len5p


def test_v3_roundtrip_preserves_loopouts(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    s5p, s3p = _two_fwd_strands(part)
    fwd2 = part.getStrandSets(2)[0]
    s_extra = fwd2.createStrand(0, 3)
    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 3)
    part.createLoopout(s3p, s3p.idx3Prime(), s_extra, s_extra.idx5Prime(), 7)
    oligo = s5p.oligo()
    total = oligo.length()
    # partially known sequence: cover only the first part, rest stays unknown
    oligo.applySequence(('ACGT' * 4)[:16])

    obj = json.loads(cnencode.encode(doc))
    # loopouts are serialized with their length; plain crossovers stay 6-tuples
    xovers = obj['parts'][0]['xovers']
    assert sorted(len(x) for x in xovers) == [7, 7]

    from cadnano.document import Document
    doc2 = Document()
    v3decode.decode(doc2, obj, emit_signals=False)
    part2 = doc2.activePart()
    assert not part2.hasLoopouts() is False
    a = part2.getStrandSets(0)[0].getStrand(0)
    b = part2.getStrandSets(1)[0].getStrand(20)
    c = part2.getStrandSets(2)[0].getStrand(0)
    assert a.connection3p() is b and a.loopoutLength() == 3
    assert b.connection3p() is c and b.loopoutLength() == 7
    assert a.oligo() is b.oligo() is c.oligo()
    assert a.oligo().length() == total
    # the known 5' prefix survives the round trip, unknowns remain unknown
    assert a.oligo().sequence()[:16] == ('ACGT' * 4)[:16]
    csv_text = part2.getSequences()
    assert '?' in csv_text


def test_legacy_export_rejects_loopouts(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    # without loopouts the new guard never fires (this minimal part fails
    # legacy export later for unrelated grid-type reasons, not loopouts)
    with pytest.raises(ValueError) as excinfo:
        cnencode.encode(doc, legacy=True)
    assert not isinstance(excinfo.value, LoopoutLegacyExportError)
    s5p, s3p = _two_fwd_strands(part)
    part.createLoopout(s5p, s5p.idx3Prime(), s3p, s3p.idx5Prime(), 2)
    with pytest.raises(LoopoutLegacyExportError):
        cnencode.encode(doc, legacy=True)


def test_plain_crossovers_keep_legacy_tuple_and_loop_ability(cnapp):
    doc = cnapp.document
    part = create3Helix(doc, [0, 0, 1], 42)
    fwd0, _ = part.getStrandSets(0)
    s = fwd0.createStrand(5, 20)
    # legacy in-strand crossover circularization still allowed for xovers
    part.createXover(s, s.idx3Prime(), s, s.idx5Prime())
    assert s.oligo().isCircular()
    xover_list = []
    fwd0.dump(xover_list)
    assert all(len(x) == 6 for x in xover_list)
    assert not part.hasLoopouts()
