import pathlib

import regex

# define some regex stuff
blanks = r"\s+"  # matches >=1  subsequent whitespace characters
sign = r"[-+]?"  # contains either a '-' or a '+' symbol or none of both
# matches ints or floats of forms '1.105' or '1.105e-5' or '1.105e5' or '1.105e+5'
float_re = r"\d+(?:\.\d+)?(?:[e][-+]?\d+)?"

tree_id_re = r"\d+"  # tree ID is an int
llh_re = rf"{sign}{float_re}"  # likelihood is a signed floating point
deltaL_re = rf"{sign}{float_re}"  # deltaL is a signed floating point
# test result entry is of form '0.123 +'
test_result_re = rf"{float_re}{blanks}{sign}"

stat_test_name = r"[a-zA-Z-]+"

# table header is of form:
# Tree      logL    deltaL  bp-RELL    p-KH     p-SH    p-WKH    p-WSH       c-ELW       p-AU
table_header = rf"Tree{blanks}logL{blanks}deltaL{blanks}(?:({stat_test_name})\s*)*"
table_header_re = regex.compile(table_header)

# a table entry in the .iqtree file looks for example like this:
# 5 -5708.931281 1.7785e-06  0.0051 -  0.498 +  0.987 +  0.498 +  0.987 +      0.05 +    0.453 +
table_entry = rf"({tree_id_re}){blanks}({llh_re}){blanks}({deltaL_re}){blanks}(?:({test_result_re})\s*)*"
table_entry_re = regex.compile(table_entry)

# if there is only a single plausible tree
# the line will look like this:
# 1  -88.9544627       0
table_entry_single_plausible_tree = (
    rf"({tree_id_re}){blanks}({llh_re}){blanks}({deltaL_re})\s*"
)
table_entry_single_plausible_tree_re = regex.compile(table_entry_single_plausible_tree)

START_STRING = "USER TREES"
END_STRING = "TIME STAMP"

TEST_NAMES = ["bp-RELL", "p-KH", "p-SH", "p-WKH", "p-WSH", "c-ELW", "p-AU"]


def get_relevant_section(input_file: pathlib.Path) -> list[str]:
    """
    Returns the section between the START_STRING and END_STRING in the given file.

    Args:
        input_file (pathlib.Path): Path to the .iqtree file to extract the section from.

    Returns:
        A list of strings, each representing a line of the relevant section.

    Raises:
        ValueError: If the section between START_STRING and END_STRING is empty.
    """
    content = input_file.read_text().splitlines()

    # now let's find the relevant lines
    # the relevant lines are only between the start and end string
    start = 0
    end = 0

    for i, line in enumerate(content):
        if START_STRING in line:
            start = i
        if END_STRING in line:
            end = i

    if start == end:
        raise ValueError(
            f"The section between START_STRING {START_STRING} and END_STRING {END_STRING} is empty. Please check the input file {input_file}."
        )

    return content[start:end]


def _get_default_entry() -> dict:
    """
    Returns a default entry for a single plausible tree.

    Returns:
        A dict containing the default entry for a single plausible tree.
    """
    return {
        "plausible": 1,
        "tests": {
            "bp-RELL": {"score": 1, "significant": True},
            "p-KH": {"score": 1, "significant": True},
            "p-SH": {"score": 1, "significant": True},
            "p-WKH": {"score": 1, "significant": True},
            "p-WSH": {"score": 1, "significant": True},
            "c-ELW": {"score": 1, "significant": True},
            "p-AU": {"score": 1, "significant": True},
        },
    }


def _regex_group_to_test_results(raw_results: list[str]) -> dict:
    assert len(TEST_NAMES) == len(raw_results)

    data = {"tests": {}}
    num_passed = 0

    for i, test in enumerate(TEST_NAMES):
        test_result = raw_results[i]
        score, significant = test_result.split(" ")
        score = score.strip()
        significant = significant.strip()
        data["tests"][test] = {}
        data["tests"][test]["score"] = float(score)
        data["tests"][test]["significant"] = True if significant == "+" else False

        if data["tests"][test]["significant"]:
            num_passed += 1

    data["plausible"] = num_passed == len(data["tests"].keys())
    return data


def get_cleaned_table_entries(table_section: list[str]) -> list[dict]:
    """
    Returns a list of dicts, each dict contains the iqtree test results for the respective tree.

    Args:
        table_section (list[str]): A list of strings, each representing a line of the relevant section.

    Returns:
        A list of dicts. Each dict contains the tree_id, llh, deltaL and all results of the performed
            iqtree tests.

    """
    entries = []
    for line in table_section:
        line = line.strip()
        # match the line against the regex defined above for a table entry
        m = regex.match(table_entry_re, line)

        # and match the line against the regex for a table entry in case of a single plausible tree
        m_single_tree = regex.match(table_entry_single_plausible_tree_re, line)

        if m:
            # transform the raw results to a python dict
            entry = _regex_group_to_test_results(m.captures(4))
            entries.append(entry)
        elif m_single_tree:
            # if a match for a truncated table entry was found: we only have a single plausible tree
            # => add the entry manually
            entry = _get_default_entry()
            entries.append(entry)
        elif "= tree" in line:
            # indicates that a tree is identical to one seen before
            # => duplicate the results of this tree
            _, id_of_identical_tree = line.rsplit(" ", 1)
            id_of_identical_tree = int(id_of_identical_tree)

            # IQ-Tree reports the results 1-indexed
            # => to get the correct results we need to subtract one and access the entries
            entry = entries[id_of_identical_tree - 1].copy()
            entries.append(entry)

    if not entries:
        raise ValueError(
            "No line in the given section matches the regex. Compare the regex and the given section. "
            "Maybe the format has changed."
        )

    return entries


def get_iqtree_results(iqtree_file: pathlib.Path) -> list[dict]:
    """
    Returns a list of dicts, each dict contains the iqtree test results for the respective tree.

    Args:
        iqtree_file (pathlib.Path): Path to the .iqtree file to extract the results from.

    Returns:
        A list of dicts. Each dict contains the tree_id, llh, deltaL and all results of the performed
            iqtree tests.
    """
    section = get_relevant_section(iqtree_file)
    entries = get_cleaned_table_entries(section)
    return entries
