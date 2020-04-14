import numpy as np
import sage_analysis.example_plots
import sage_analysis.example_calcs
from sage_analysis.utils import generate_func_dict, select_random_indices
import types

def test_func_dict_calcs():

    plot_toggles = {"SMF": 1}
    module_name = "sage_analysis.example_calcs"
    function_prefix = "calc_"
    func_dict = generate_func_dict(plot_toggles, module_name, function_prefix)

    assert list(func_dict.keys()) == ["SMF"]

    assert isinstance(func_dict["SMF"], tuple)
    assert isinstance(func_dict["SMF"][0], types.FunctionType)
    assert isinstance(func_dict["SMF"][1], dict)
    assert len(func_dict["SMF"][1]) == 0

def test_func_dict_plots():

    plot_toggles = {"SMF": 1, "quiescent": 1}
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

# Use pytest mark
def test_random_indices_less_than_number_available():
    """
    Request less than the number of inds available across all files, but more than is in this file -> returns a random
    subset of the indices.
    """

    seed = 666
    inds = np.arange(10)
    global_num_inds_available = 100
    global_num_inds_requested = 50
    inds = list(select_random_indices(inds, global_num_inds_available, global_num_inds_requested, seed))
    assert inds == [2, 6, 9, 4, 3]
