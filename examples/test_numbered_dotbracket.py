from numbered_dotbracket import NumberedDotbracket

def test_without_substr():
    db = "([)(([)])(([)({)])([)](([)](]))"
    a = NumberedDotbracket(db, list(range(len(db))))
    removed, remaining = a.without_substr("([)]")
    assert remaining == "([)()(([)({)])((]))"
    assert removed[0] == NumberedDotbracket("([)]", [4,5,6,7])
    assert len(removed)==3

def test_without_short_helices():
    db = "(((.((...)))))"
    a = NumberedDotbracket(db, list(range(len(db))))
    a2 = a.without_short_helices(3)
    assert a2 == "(((........)))"
    a3 = a.without_short_helices(4)
    assert a3 == ".............."
