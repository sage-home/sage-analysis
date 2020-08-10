from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

from tqdm import tqdm

from sage_analysis.model import Model


class DataClass(ABC):
    """
    An abstract baseclass for handling the various **SAGE** output formats. It should not be instantiated directly;
    instead the underlying subclasses for each format (e.g., :py:class:`~sage_analysis.sage_binary.SageBinaryData` and
    :py:class:`~sage_analysis.sage_binary.SageBinaryData`).

    We refer to :doc:`../user/data_class` for more information about adding your own Data Class to ingest custom data
    formats.
    """

    @abstractmethod
    def determine_volume_analyzed(self, model: Model, **kwargs: Any) -> float:
        """
        Determines the volume analyzed. This can be smaller than the total simulation box.

        Parameters
        ----------
        model : :py:class:`~sage_analysis.model.Model` instance
            The model that this data class is associated with.

        **kwargs : any
            Extra arguments to allow other data classes to pass extra arguments to their implementation.

        Returns
        -------
        volume : float
            The numeric volume being processed during this run of the code in (Mpc/h)^3.
        """
        pass  # pragma: no cover

    @abstractmethod
    def read_sage_params(self, sage_file_path: str, *kwargs: Any) -> Dict[str, Any]:
        """
        Read the **SAGE** parameter file.

        Parameters
        ----------
        sage_file_path: string
            Path to the **SAGE** parameter file.

        **kwargs : any
            Extra arguments to allow other data classes to pass extra arguments to their implementation.

        Returns
        -------
        model_dict: dict [str, var]
            Dictionary containing the parameter names and their values.
        """
        pass  # pragma: no cover

    @abstractmethod
    def determine_num_gals(self, model: Model, **kwargs: Any):
        """
        Determines the number of galaxies in all files for this :py:class:`~sage_analysis.model.Model`.

        Parameters
        ----------
        model: :py:class:`~sage_analysis.model.Model` class
            The :py:class:`~sage_analysis.model.Model` we're reading data for.

        **kwargs : any
            Extra arguments to allow other data classes to pass extra arguments to their implementation.
        """
        pass  # pragma: no cover

    @abstractmethod
    def read_gals(
        self,
        model: Model,
        file_num: int,
        snapshot: int,
        pbar: Optional[tqdm] = None,
        plot_galaxies: bool = False,
        debug: bool = False,
        **kwargs: Any,
    ) -> Any:
        """
        Reads the galaxies of a model file for the specified file number and snapshot.

        Parameters
        ----------
        model : :py:class:`~sage_analysis.model.Model` class
            The :py:class:`~sage_analysis.model.Model` we're reading data for.

        file_num : int
            Suffix number of the file we're reading.

        snapshot : int
            The snapshot we're reading.

        pbar : ``tqdm`` class instance, optional
            Bar showing the progress of galaxy reading.  If not specified, progress bar will not show.

        plot_galaxies : bool, optional
            If specified, plots and saves the 3D distribution of galaxies for this file.

        debug : bool, optional
            If specified, prints out extra useful debug information.

        **kwargs : any
            Extra arguments to allow other data classes to pass extra arguments to their implementation.

        Returns
        -------
        gals : The format is specified by the underlying data class implementation
            The galaxies for this file.
        """
        pass  # pragma: no cover

    @abstractmethod
    def update_snapshot_and_data_path(self, model: Model, snapshot: int, **kwargs: Any) -> None:
        """
        Updates the :py:attr:`~sage_analysis.model.Model._sage_data_path` to point to a new redshift file (if necessary
        for the underlying implementation). Uses the redshift array :py:attr:`~sage_analysis.model.Model.redshifts`.

        Parameters
        ----------
        snapshot : int
            Snapshot we're updating :py:attr:`~sage_analysis.model.Model._sage_data_path` to
            point to.

        **kwargs : any
            Extra arguments to allow other data classes to pass extra arguments to their implementation.
        """
        pass  # pragma: no cover

    @abstractmethod
    def close_file(self, model: Model, **kwargs) -> None:
        """
        Closes an open galaxy file.  This is useful when reading the HDF5 data format where a single file contains many
        snapshots. For the binary format, this is an empty method.

        Parameters
        ----------
        **kwargs : any
            Extra arguments to allow other data classes to pass extra arguments to their implementation.
        """
        pass  # pragma: no cover
