from labelgenerator.iqtree_parser import get_iqtree_results


def test_get_iqtree_results_dna_msa(log_dir):
    expected_result = [
        {
            "tests": {
                "bp-RELL": {"score": 0.0012, "significant": False},
                "p-KH": {"score": 0.466, "significant": True},
                "p-SH": {"score": 0.945, "significant": True},
                "p-WKH": {"score": 0.466, "significant": True},
                "p-WSH": {"score": 0.993, "significant": True},
                "c-ELW": {"score": 0.0431, "significant": True},
                "p-AU": {"score": 0.787, "significant": True},
            },
            "plausible": False,
        },  # 1
        {
            "tests": {
                "bp-RELL": {"score": 0.0046, "significant": False},
                "p-KH": {"score": 0.615, "significant": True},
                "p-SH": {"score": 1.0, "significant": True},
                "p-WKH": {"score": 0.534, "significant": True},
                "p-WSH": {"score": 0.984, "significant": True},
                "c-ELW": {"score": 0.0452, "significant": True},
                "p-AU": {"score": 0.626, "significant": True},
            },
            "plausible": False,
        },  # 2
        {
            "tests": {
                "bp-RELL": {"score": 0.124, "significant": True},
                "p-KH": {"score": 0.232, "significant": True},
                "p-SH": {"score": 0.541, "significant": True},
                "p-WKH": {"score": 0.239, "significant": True},
                "p-WSH": {"score": 0.41, "significant": True},
                "c-ELW": {"score": 0.0446, "significant": True},
                "p-AU": {"score": 0.17, "significant": True},
            },
            "plausible": True,
        },  # 3
        {
            "tests": {
                "bp-RELL": {"score": 0.0082, "significant": False},
                "p-KH": {"score": 0.464, "significant": True},
                "p-SH": {"score": 0.945, "significant": True},
                "p-WKH": {"score": 0.0453, "significant": False},
                "p-WSH": {"score": 0.179, "significant": True},
                "c-ELW": {"score": 0.043, "significant": True},
                "p-AU": {"score": 0.786, "significant": True},
            },
            "plausible": False,
        },  # 4
        {
            "tests": {
                "bp-RELL": {"score": 0.0125, "significant": False},
                "p-KH": {"score": 0.464, "significant": True},
                "p-SH": {"score": 0.945, "significant": True},
                "p-WKH": {"score": 0.0001, "significant": False},
                "p-WSH": {"score": 0.0001, "significant": False},
                "c-ELW": {"score": 0.043, "significant": True},
                "p-AU": {"score": 0.788, "significant": True},
            },
            "plausible": False,
        },  # 5
        {
            "tests": {
                "bp-RELL": {"score": 0.0186, "significant": True},
                "p-KH": {"score": 0.464, "significant": True},
                "p-SH": {"score": 0.945, "significant": True},
                "p-WKH": {"score": 0.0479, "significant": False},
                "p-WSH": {"score": 0.181, "significant": True},
                "c-ELW": {"score": 0.043, "significant": True},
                "p-AU": {"score": 0.787, "significant": True},
            },
            "plausible": False,
        },  # 6
        {
            "tests": {
                "bp-RELL": {"score": 0.028, "significant": True},
                "p-KH": {"score": 0.458, "significant": True},
                "p-SH": {"score": 0.943, "significant": True},
                "p-WKH": {"score": 0.201, "significant": True},
                "p-WSH": {"score": 0.725, "significant": True},
                "c-ELW": {"score": 0.0423, "significant": True},
                "p-AU": {"score": 0.783, "significant": True},
            },
            "plausible": True,
        },  # 7
        {
            "tests": {
                "bp-RELL": {"score": 0.036, "significant": True},
                "p-KH": {"score": 0.464, "significant": True},
                "p-SH": {"score": 0.945, "significant": True},
                "p-WKH": {"score": 0.0114, "significant": False},
                "p-WSH": {"score": 0.047, "significant": False},
                "c-ELW": {"score": 0.0429, "significant": True},
                "p-AU": {"score": 0.787, "significant": True},
            },
            "plausible": False,
        },  # 8
        {
            "tests": {
                "bp-RELL": {"score": 0.0309, "significant": True},
                "p-KH": {"score": 0.453, "significant": True},
                "p-SH": {"score": 0.846, "significant": True},
                "p-WKH": {"score": 0.453, "significant": True},
                "p-WSH": {"score": 0.963, "significant": True},
                "c-ELW": {"score": 0.0731, "significant": True},
                "p-AU": {"score": 0.553, "significant": True},
            },
            "plausible": True,
        },  # 9
        {
            "tests": {
                "bp-RELL": {"score": 0.0582, "significant": True},
                "p-KH": {"score": 0.452, "significant": True},
                "p-SH": {"score": 0.845, "significant": True},
                "p-WKH": {"score": 0.254, "significant": True},
                "p-WSH": {"score": 0.726, "significant": True},
                "c-ELW": {"score": 0.0729, "significant": True},
                "p-AU": {"score": 0.558, "significant": True},
            },
            "plausible": True,
        },  # 10
        {
            "tests": {
                "bp-RELL": {"score": 0.027, "significant": True},
                "p-KH": {"score": 0.424, "significant": True},
                "p-SH": {"score": 0.961, "significant": True},
                "p-WKH": {"score": 0.424, "significant": True},
                "p-WSH": {"score": 0.952, "significant": True},
                "c-ELW": {"score": 0.0452, "significant": True},
                "p-AU": {"score": 0.629, "significant": True},
            },
            "plausible": True,
        },  # 11
        {
            "tests": {
                "bp-RELL": {"score": 0.074, "significant": True},
                "p-KH": {"score": 0.453, "significant": True},
                "p-SH": {"score": 0.846, "significant": True},
                "p-WKH": {"score": 0.453, "significant": True},
                "p-WSH": {"score": 0.976, "significant": True},
                "c-ELW": {"score": 0.0731, "significant": True},
                "p-AU": {"score": 0.557, "significant": True},
            },
            "plausible": True,
        },  # 12
        {
            "tests": {
                "bp-RELL": {"score": 0.165, "significant": True},
                "p-KH": {"score": 0.437, "significant": True},
                "p-SH": {"score": 0.831, "significant": True},
                "p-WKH": {"score": 0.437, "significant": True},
                "p-WSH": {"score": 0.844, "significant": True},
                "c-ELW": {"score": 0.0822, "significant": True},
                "p-AU": {"score": 0.438, "significant": True},
            },
            "plausible": True,
        },  # 13
        {
            "tests": {
                "bp-RELL": {"score": 0.106, "significant": True},
                "p-KH": {"score": 0.452, "significant": True},
                "p-SH": {"score": 0.846, "significant": True},
                "p-WKH": {"score": 0.346, "significant": True},
                "p-WSH": {"score": 0.98, "significant": True},
                "c-ELW": {"score": 0.073, "significant": True},
                "p-AU": {"score": 0.56, "significant": True},
            },
            "plausible": True,
        },  # 14
        {
            "tests": {
                "bp-RELL": {"score": 0.0418, "significant": True},
                "p-KH": {"score": 0.385, "significant": True},
                "p-SH": {"score": 0.979, "significant": True},
                "p-WKH": {"score": 0.385, "significant": True},
                "p-WSH": {"score": 0.967, "significant": True},
                "c-ELW": {"score": 0.0452, "significant": True},
                "p-AU": {"score": 0.626, "significant": True},
            },
            "plausible": True,
        },  # 15
        {
            "tests": {
                "bp-RELL": {"score": 0.0183, "significant": False},
                "p-KH": {"score": 0.0823, "significant": True},
                "p-SH": {"score": 0.14, "significant": True},
                "p-WKH": {"score": 0.0673, "significant": True},
                "p-WSH": {"score": 0.113, "significant": True},
                "c-ELW": {"score": 0.0101, "significant": False},
                "p-AU": {"score": 0.0201, "significant": False},
            },
            "plausible": False,
        },  # 16
        {
            "tests": {
                "bp-RELL": {"score": 0.0513, "significant": True},
                "p-KH": {"score": 0.179, "significant": True},
                "p-SH": {"score": 0.927, "significant": True},
                "p-WKH": {"score": 0.015, "significant": False},
                "p-WSH": {"score": 0.108, "significant": True},
                "c-ELW": {"score": 0.045, "significant": True},
                "p-AU": {"score": 0.632, "significant": True},
            },
            "plausible": False,
        },  # 17
        {
            "tests": {
                "bp-RELL": {"score": 0.0618, "significant": True},
                "p-KH": {"score": 0.022, "significant": False},
                "p-SH": {"score": 0.927, "significant": True},
                "p-WKH": {"score": 0.0, "significant": False},
                "p-WSH": {"score": 0.0074, "significant": False},
                "c-ELW": {"score": 0.045, "significant": True},
                "p-AU": {"score": 0.628, "significant": True},
            },
            "plausible": False,
        },  # 18
        {
            "tests": {
                "bp-RELL": {"score": 0.0705, "significant": True},
                "p-KH": {"score": 0.0862, "significant": True},
                "p-SH": {"score": 0.944, "significant": True},
                "p-WKH": {"score": 0.0862, "significant": True},
                "p-WSH": {"score": 0.395, "significant": True},
                "c-ELW": {"score": 0.0451, "significant": True},
                "p-AU": {"score": 0.625, "significant": True},
            },
            "plausible": True,
        },  # 19
        {
            "tests": {
                "bp-RELL": {"score": 0.0626, "significant": True},
                "p-KH": {"score": 0.465, "significant": True},
                "p-SH": {"score": 0.945, "significant": True},
                "p-WKH": {"score": 0.0788, "significant": True},
                "p-WSH": {"score": 0.291, "significant": True},
                "c-ELW": {"score": 0.0431, "significant": True},
                "p-AU": {"score": 0.787, "significant": True},
            },
            "plausible": True,
        },  # 20
    ]

    result = get_iqtree_results(log_dir / "DNA.iqtree.iqtree")
    assert result == expected_result


def test_get_iqtree_results_aa_msa(log_dir):
    expected_result = [
        {
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
    ] * 20

    result = get_iqtree_results(log_dir / "AA.iqtree.iqtree")
    assert result == expected_result


def test_get_iqtree_results_morph_msa(log_dir):
    entry_tree_1 = {
        "plausible": 1,
        "tests": {
            "bp-RELL": {"score": 0.505, "significant": True},
            "p-KH": {"score": 0.902, "significant": True},
            "p-SH": {"score": 1, "significant": True},
            "p-WKH": {"score": 0.902, "significant": True},
            "p-WSH": {"score": 0.902, "significant": True},
            "c-ELW": {"score": 0.5, "significant": True},
            "p-AU": {"score": 0.84, "significant": True},
        },
    }

    entry_tree_4 = {
        "plausible": 1,
        "tests": {
            "bp-RELL": {"score": 0.495, "significant": True},
            "p-KH": {"score": 0.0983, "significant": True},
            "p-SH": {"score": 0.0983, "significant": True},
            "p-WKH": {"score": 0.0983, "significant": True},
            "p-WSH": {"score": 0.0983, "significant": True},
            "c-ELW": {"score": 0.5, "significant": True},
            "p-AU": {"score": 0.16, "significant": True},
        },
    }

    expected_result = [
        entry_tree_1,  # 1
        entry_tree_1,  # 2
        entry_tree_1,  # 3
        entry_tree_4,  # 4
        entry_tree_4,  # 5
        entry_tree_4,  # 6
        entry_tree_1,  # 7
        entry_tree_1,  # 8
        entry_tree_4,  # 9
        entry_tree_4,  # 10
        entry_tree_4,  # 11
        entry_tree_1,  # 12
        entry_tree_1,  # 13
        entry_tree_4,  # 14
        entry_tree_4,  # 15
        entry_tree_1,  # 16
        entry_tree_4,  # 17
        entry_tree_4,  # 18
        entry_tree_4,  # 19
        entry_tree_1,  # 20
    ]

    result = get_iqtree_results(log_dir / "MORPH.iqtree.iqtree")
    assert result == expected_result
