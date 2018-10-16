from string import ascii_letters
from copy import copy

from hypothesis.strategies import composite, text, sampled_from, one_of, data
from hypothesis.strategies import characters, integers, none, just, lists, booleans
from hypothesis import given, assume
import forgi.graph.sequence as fgs
import forgi.graph.residue as fgr

@composite
def get_insertion_code(draw, min_codepoint=None):
    if min_codepoint == None:
        min_codepoint=65
    return draw(characters(min_codepoint=min_codepoint, max_codepoint=90))

@composite
def resid_strategy(draw, i_kwargs={}):
    """
    strategy for generating a resid.

    Look at the code, to see what we currently support as resid.
    """
    chain = draw(one_of(none(), text(min_size=1, alphabet=ascii_letters)))
    het_flag=" " # Currently, a non-empty het-flag is not supported!
    resid = draw(integers(**i_kwargs).filter(lambda x: x!=0))
    insertion_code = draw(one_of(just(" "), get_insertion_code()))
    return fgr.RESID(chain, (het_flag, resid, insertion_code))

@composite
def next_residue(draw, last_residue, with_missing=False, with_insertion=False, with_chainbreak=False):
    assume (with_missing+with_insertion+with_chainbreak<=1)
    if with_chainbreak:
        assume (last_residue[0] is not None) # Don't sample chain-breaks without chain
        next_res = draw(resid_strategy(
                            i_kwargs={"min_value":-25,"max_value":200}).filter(
                            lambda x: x[0] is not None and x[0]>last_residue[0]))
        return next_res

    next_i = last_residue[1][1]
    if not with_insertion:
        next_i+=1
        i_code = " "
    else:
        prev_insertion = last_residue[1][2]
        if prev_insertion !=" ":
            mcp=ord(prev_insertion)+1
        else:
            mcp=None
        assume (mcp is None or mcp<90)
        i_code = draw(get_insertion_code(min_codepoint=mcp))
    if with_missing:
        next_i+=draw(integers(min_value=1, max_value=400))
    return fgr.RESID(last_residue[0], (" ", next_i, i_code))

@composite
def seqids_strategy(draw, with_missing=False, with_insertions=False, with_chainbreak=False):
    if not with_missing and not with_insertions and not with_chainbreak:
        return draw(one_of(just([]), _seqids_strategy_core))
    else:
        return draw(_seqids_strategy_core(with_missing, with_insertions, with_chainbreak))

@composite
def _seqids_strategy_core(draw, with_missing, with_insertions, with_chainbreak):
    first_res = draw(resid_strategy(i_kwargs={"min_value":-25, "max_value":200}))
    if with_chainbreak:
        assume (first_res[0] is not None)
    length = draw(integers(min_value=with_missing+with_insertions+with_chainbreak+1, max_value=100))
    # Force with_missing, with_chainbreak and with_insertions to get triggered at least once!
    if with_missing:
        missing_index = draw(integers(min_value=1, max_value=length-1))
    else:
        missing_index=-1
    if with_insertions:
        insertion_index = draw(integers(min_value=1, max_value=length-1).filter(
                                lambda x:x!=missing_index))
    else:
        insertion_index = -1
    if with_chainbreak:
        chainbreak_index = draw(integers(min_value=1, max_value=length-1).filter(
                                lambda x:x!=missing_index and x!=insertion_index))
    else:
        insertion_index = -1
    seq_ids = [first_res]
    for i in range(1, length):
        if i==missing_index:
            wm=True
        else:
            wm=False
        if i==insertion_index:
            wi=True
        else:
            wi=False
        if i==chainbreak_index:
            wc=True
        else:
            wc=False
        extra_true = draw(integers(max_value=10))>8
        if extra_true and wm+wc+wi==0:
            choice = draw(sampled_from([0,1,2]))
            if choice==0:
                wm=True
            elif choice==1:
                wi=True
            else:
                wc=True
        seq_ids.append(draw(next_residue(seq_ids[-1], wm, wi, wc)))
    return seq_ids

@composite
def lists_with_duplicates(draw, elem_strategy):
    l = draw(lists(elem_strategy, min_size=1))
    elem = draw(sampled_from(l))
    i = draw(integers(min_value=0, max_value=len(l)-1))
    l.insert(i, copy(elem))
    return l

@given(resid=resid_strategy(), resname=sampled_from("AUGC"))
def test_missing_residue(resid, resname):
    #print("testing", resid, resname)
    mr = fgs.MissingResidue(resid, resname)
    mr_roundtrip = fgs.MissingResidue.from_bg_fields(mr.to_bg_string().split())
    assert mr.resid == mr_roundtrip.resid
    assert mr.res_name == mr_roundtrip.res_name

@given(seq=lists(resid_strategy(), unique=True))
def test_seqidList_index_and_getitem_consistent(seq):
    seqlist = fgs.SeqidList(seq)
    for elem in seqlist:
        assert seqlist[seqlist.index(elem)]==elem

@given(seq1=lists(resid_strategy(), unique=True), seq2=lists(resid_strategy(), unique=True))
def test_seqlist_equal(seq1, seq2):
    assume (seq1 != seq2)
    assert fgs.SeqidList(seq1) != fgs.SeqidList(seq2)
    assert fgs.SeqidList(seq1) == fgs.SeqidList(seq1)
    assert fgs.SeqidList(seq2) == fgs.SeqidList(seq2)

@given(lists_with_duplicates(resid_strategy()))
def test_seqids_duplicate(seq):
    try:
        fgs.SeqidList(seq)
    except ValueError:
        pass
    else:
        assert False, "ValueError not raised for duplicate seqid {} in {}".format(elem, seq)

@given(seqids_strategy(True, True, True))
def test_seqids_strategy(seqids):
    assert len(seqids)>=3
    found_c = False
    found_m = False
    found_i = False
    for i, s in enumerate(seqids):
        if i>0:
            assert seqids[i-1]<s
            if seqids[i-1][1][1]+1<s[1][1]:
                found_m = True
            if s[1][2]!=" ":
                found_i=True
                if seqids[i-1][0]==s[0]:
                    assert seqids[i-1][1][1]==s[1][1], "{}, {}".format(seqids[i-1], s)
            if s[0]!=seqids[i-1][0]:
                found_c=True
    assert found_c
    assert found_i
    assert found_m
