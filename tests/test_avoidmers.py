import pytest
import json
from pathlib import Path
import csv

parent = Path(__file__).resolve().parent
answers_path = parent.joinpath("answers.json")

# with open(answers_path, 'r') as f:
#     data = json.load(f)

data = {
  "test4": [
            ("agatag", 0, 6),
            ("gatagat", 1, 8),
            ("atagat", 2, 8),
            ("tagatatg", 3, 11),
            ("agatatg", 4, 11),
            ("gatatg", 5, 11),
        ],
  "test2": [
            ("catgcagtgacgatgctagcaacaatcatgactta", 0, 35)
            ],
  "test3": [
            ("agagtcagttgtatag", 0, 16),
            ("gagtcagttgtataga", 1, 17),
  ]
}

@pytest.mark.parametrize("test_input", ["test3", "test4"])
def test_eval(test_input):
    global data

    test_suite = data[test_input]
    results = []
    path = parent.parent.joinpath("avoidmers",
                                     "pattern_extractions_abacaba")

    if test_input == "test4":
        file = path.joinpath(f"test_abacaba_words_length_6_seq_{test_input}.txt")
    else:
        file = path.joinpath(f"test_abacaba_words_length_6_seq_{test_input}.maximal.txt")

    with open(file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            results.append((row["sequence"], int(row["start"]), int(row["end"])))

    assert len(results) == len(test_suite)
    assert set(results) == set(test_suite)
