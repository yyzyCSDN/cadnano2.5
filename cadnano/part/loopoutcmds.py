# -*- coding: utf-8 -*-
"""Undo commands for single-stranded loopouts.

A loopout is a connection from the 3' end of ``strand5p`` to the 5' end of
``strand3p`` carrying one or more unpaired bases that occupy no virtual helix
grid point. A connection with loopout length 0 is a plain crossover.
"""
from cadnano.proxies.cnproxy import UndoCommand
from cadnano.strand import Strand
from cadnano.cntypes import (
    NucleicAcidPartT,
    StrandT
)


class CreateLoopoutCommand(UndoCommand):
    """Create a loopout from the 3' end of ``strand5p`` to the 5' end of
    ``strand3p`` with ``loopout_length`` unpaired bases.

    Mirrors :class:`CreateXoverCommand` oligo bookkeeping, additionally
    accounting for the unpaired loopout bases. The bound strands, their grid
    coordinates and existing insertions/skips are never modified.
    """

    def __init__(self,  part: NucleicAcidPartT,
                        strand5p: StrandT, strand5p_idx: int,
                        strand3p: StrandT, strand3p_idx: int,
                        loopout_length: int,
                        update_oligo: bool = True):
        super(CreateLoopoutCommand, self).__init__("create loopout")
        self._part = part
        self._strand5p = strand5p
        self._strand5p_idx = strand5p_idx
        self._strand3p = strand3p
        self._strand3p_idx = strand3p_idx
        self._loopout_length = loopout_length
        self._old_oligo3p = strand3p.oligo()
        self._update_oligo = update_oligo
    # end def

    def redo(self):
        strand5p = self._strand5p
        strand3p = self._strand3p
        olg5p = strand5p.oligo()
        old_olg3p = self._old_oligo3p

        doc = strand5p.document()
        doc.removeStrandFromSelection(strand5p)
        doc.removeStrandFromSelection(strand3p)

        fSetOligo = Strand.setOligo
        if self._update_oligo:
            if olg5p == strand3p.oligo():
                olg5p._setLoop(True)
            else:
                # the merged oligo gains the 3' oligo plus the unpaired bases
                olg5p._incrementLength(old_olg3p.length() + self._loopout_length,
                                        emit_signals=True)
                old_olg3p.removeFromPart(emit_signals=True)
                for strand in strand3p.generator3pStrand():
                    fSetOligo(strand, olg5p, emit_signals=True)

        # install the connection and the single-stranded loopout
        strand5p.setConnection3p(strand3p)
        strand3p.setConnection5p(strand5p)
        strand5p.setLoopout(self._loopout_length)

        if self._update_oligo:
            strand5p.strandConnectionChangedSignal.emit(strand5p)
            strand3p.strandConnectionChangedSignal.emit(strand3p)
    # end def

    def undo(self):
        part = self._part
        strand5p = self._strand5p
        strand3p = self._strand3p
        old_olg3p = self._old_oligo3p
        olg5p = strand5p.oligo()

        doc = strand5p.document()
        doc.removeStrandFromSelection(strand5p)
        doc.removeStrandFromSelection(strand3p)

        # uninstall the loopout and its connection
        strand5p.setLoopout(0)
        strand5p.setConnection3p(None)
        strand3p.setConnection5p(None)

        fSetOligo = Strand.setOligo
        if self._update_oligo:
            if old_olg3p.isCircular():
                old_olg3p._setLoop(False)
            else:
                olg5p._decrementLength(old_olg3p.length() + self._loopout_length,
                                       emit_signals=True)
                old_olg3p.addToPart(part, emit_signals=True)
                for strand in strand3p.generator3pStrand():
                    fSetOligo(strand, old_olg3p, emit_signals=True)

            strand5p.strandConnectionChangedSignal.emit(strand5p)
            strand3p.strandConnectionChangedSignal.emit(strand3p)
            part.reemitActiveVirtualHelix()
    # end def
# end class


class ChangeLoopoutLengthCommand(UndoCommand):
    """Change the length of an existing loopout connection (length 0 being a
    plain crossover), preserving assigned bases on both bound strands and the
    5' prefix of the loopout sequence.
    """

    def __init__(self,  strand5p: StrandT,
                        new_length: int,
                        update_oligo: bool = True):
        super(ChangeLoopoutLengthCommand, self).__init__("change loopout length")
        assert strand5p.connection3p() is not None
        self._strand5p = strand5p
        self._strand3p = strand5p.connection3p()
        self._old_length = strand5p.loopoutLength()
        self._new_length = new_length
        self._old_sequence = strand5p._loopout_sequence
        self._update_oligo = update_oligo
    # end def

    def redo(self):
        strand5p = self._strand5p
        strand3p = self._strand3p
        strand5p.resizeLoopout(self._new_length)
        if self._update_oligo:
            strand5p.oligo()._incrementLength(self._new_length - self._old_length,
                                              emit_signals=True)
            strand5p.strandConnectionChangedSignal.emit(strand5p)
            strand3p.strandConnectionChangedSignal.emit(strand3p)
    # end def

    def undo(self):
        strand5p = self._strand5p
        strand3p = self._strand3p
        strand5p.setLoopout(self._old_length, self._old_sequence)
        if self._update_oligo:
            strand5p.oligo()._decrementLength(self._new_length - self._old_length,
                                              emit_signals=True)
            strand5p.strandConnectionChangedSignal.emit(strand5p)
            strand3p.strandConnectionChangedSignal.emit(strand3p)
    # end def
# end class
