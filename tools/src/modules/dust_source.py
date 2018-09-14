class DustSource:
    """The DustSource class is used to consider the radiation of the dust
    grains in rat calculations.
    """

    def __init__(self, file_io, parse_args):
        """Initialisation of the dust source parameters.

        Args:
            file_io : Handles file input/output and all
                necessary paths.
        """
        self.file_io = file_io
        self.parse_args = parse_args

        #: dict: Parameters of the dust source
        self.parameter = {
            'nr_photons': 0,
        }

        if self.parse_args.nr_photons_dust is not None:
            self.parameter['nr_photons'] = int(self.parse_args.nr_photons_dust)
        elif self.parse_args.simulation_type in ['rat']:
            self.parameter['nr_photons'] = 100

    def get_command_line(self):
        """Provides dust source command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the dust source.
        """
        return '\t<source_dust nr_photons = "' \
               + str(self.parameter['nr_photons']) + '">\n'

    def get_command(self):
        """Provides dust source command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider the dust source.
        """
        return self.get_command_line()
