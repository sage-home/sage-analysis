import pytest
from typing import List
import numpy as np
import sage_analysis.example_plots
import sage_analysis.example_calcs
from sage_analysis.utils import generate_func_dict, select_random_indices
import types


def test_func_dict_calcs():

    plot_toggles = {"SMF": True}
    module_name = "sage_analysis.example_calcs"
    function_prefix = "calc_"
    func_dict = generate_func_dict(plot_toggles, module_name, function_prefix)

    assert list(func_dict.keys()) == ["SMF"]

    assert isinstance(func_dict["SMF"], tuple)
    assert isinstance(func_dict["SMF"][0], types.FunctionType)
    assert isinstance(func_dict["SMF"][1], dict)
    assert len(func_dict["SMF"][1]) == 0

def test_func_dict_plots():

    plot_toggles = {"SMF": True, "quiescent": True}
    module_name = "sage_analysis.example_plots"
    function_prefix = "plot_"
    keyword_args = {
        "SMF": {"plot_sub_populations": True},
        "quiescent": {"plot_output_format": "pdf", "plot_sub_populations": True}
    }
    func_dict = generate_func_dict(plot_toggles, module_name, function_prefix, keyword_args)

    assert list(func_dict.keys()) == ["SMF", "quiescent"]

    assert isinstance(func_dict["SMF"], tuple)
    assert isinstance(func_dict["SMF"][0], types.FunctionType)
    assert isinstance(func_dict["SMF"][1], dict)
    assert func_dict["SMF"][1]["plot_sub_populations"]

    assert isinstance(func_dict["quiescent"], tuple)
    assert isinstance(func_dict["quiescent"][0], types.FunctionType)
    assert isinstance(func_dict["quiescent"][1], dict)
    assert func_dict["quiescent"][1]["plot_sub_populations"]
    assert func_dict["quiescent"][1]["plot_output_format"] == "pdf"


def test_func_dict_not_imported_module():
    """
    The functions that we try to use in the dict is in a module that hasn't been imported.
    """

    plot_toggles = {"SMF": True}
    module_name = "not_a_module.funcs"
    function_prefix = "calc_"

    with pytest.raises(KeyError):
        func_dict = generate_func_dict(plot_toggles, module_name, function_prefix)


@pytest.mark.parametrize(
    "local_inds, global_num_inds_available, global_num_inds_requested, expected_inds",
    [
        (np.arange(10), 100, 50, [2, 6, 9, 4, 3]),
        (np.arange(30), 100, 10, [12, 2, 13]),
        (np.arange(10), 100, 500, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    ]
)
def test_random_indices_less_than_number_available(
    local_inds: np.ndarray,
    global_num_inds_available: int,
    global_num_inds_requested: int,
    expected_inds: List[int],
):
    """
    Parameterize Combinations
    -------------------------
    """

    seed = 666
    inds = list(select_random_indices(local_inds, global_num_inds_available, global_num_inds_requested, seed))
    assert inds == expected_inds


