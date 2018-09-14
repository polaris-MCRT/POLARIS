#!/usr/bin/env python
# -*- coding: utf-8 -*-

import platform
from custom.server import *
from modules.base import Server


class ServerChooser:
    """The ServerChooser class provides the composition of the chosen dust
    species.
    """

    def __init__(self, parse_args, tool_type):
        """Initialisation of all usable options.

        Notes:
            To create your own dust composition, add its name to the dictionary
            and write a class with its options as a derived class of class Dust.

        Args:
            parse_args : Provides all parameters chosen
                by user when executing PolarisTools.
            tool_type (str): The toolkit that executes this function.
                'plot'  : polaris-plot
                'run'   : polaris-run
                'remote': polaris-remote
                'grid'  : polaris-gen
        """
        self.parse_args = parse_args

        self.tool_type = tool_type

        #: dict: Dictionary with all usable server/cluster
        self.server_dict = {
            'nec': NECServer,
            'rz_cluster': RZServer,
            'astro_cluster': AstroServer,
            'local_pc': Server,
        }
        update_server_dict(self.server_dict)

    def get_module(self):
        """Chooses server class from user input

        Note:
            Parameters set by PolarisTools overwrite preset values in the
            separate server classes.

        Returns:
            Instance of chosen server/cluster.
        """
        if self.tool_type == 'remote':
            if self.parse_args.server_name in self.server_dict.keys():
                server = self.server_dict[self.parse_args.server_name](self.parse_args)
            elif self.parse_args.server_name is None:
                server = self.server_dict['nec'](self.parse_args)
            else:
                raise ValueError('Server not known! You can add a new server in server.py')
            # Overwrite preset parameters from user input
            if self.parse_args.user_id is not None:
                server.parameter['user_id'] = self.parse_args.user_id
        elif self.tool_type == 'run':
            for server_name in self.server_dict.keys():
                server = self.server_dict[server_name](self.parse_args)
                if server.parameter['node_name'] != None and server.parameter['node_name'] in platform.node():
                    break
                else:
                    server = self.server_dict['local_pc'](self.parse_args)
            # Overwrite preset parameters from user input
            if self.parse_args.nr_threads is not None:
                server.parameter['nr_threads'] = self.parse_args.nr_threads
            if self.parse_args.ram_usage:
                server.parameter['ram_usage'] = self.parse_args.ram_usage
        else:
            server = self.server_dict['local_pc'](self.parse_args)
        return server


class NECServer(Server):
    """This is the server class for the NEC cluster at the Kiel university
    """

    def __init__(self, parse_args):
        """Initialisation of the server/cluster parameters.
        """
        Server.__init__(self, parse_args)

        self.parameter['node_name'] = 'nesh'
        self.parameter['address'] = 'nesh-fe.rz.uni-kiel.de:~/'
        self.parameter['queue_system'] = 'PBS'
        self.parameter['short_walltime'] = '2:00:00'
        self.parameter['short_batch_class'] = 'clexpress'
        self.parameter['medium_walltime'] = '48:00:00'
        self.parameter['medium_batch_class'] = 'clmedium'
        self.parameter['long_walltime'] = '100:00:00'
        self.parameter['long_batch_class'] = 'cllong'

    def get_command_line(self):
        """Provides server/cluster command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider a server/cluster.
        """
        new_command_line = str()
        new_command_line += '#' + self.parameter['queue_system'] + ' -o ' \
                            + self.parameter['simulation_directory'] + 'POLARIS.out' + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -j o\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -l elapstim_req=' \
                            + self.parameter['walltime'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -l memsz_job=' \
                            + self.parameter['ram_usage'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -b ' + self.parameter['node_number'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -l cpunum_job=' \
                            + str(self.parameter['nr_threads']) + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -q ' + self.parameter['batch_class'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -N ' \
                            + self.parse_args.simulation_name[0:15] + '\n'
        return new_command_line


class RZServer(Server):
    """This is the server class for the RZ-cluster at the Kiel university
    """

    def __init__(self, parse_args):
        """Initialisation of the server/cluster parameters.
        """
        Server.__init__(self, parse_args)

        self.parameter['node_name'] = 'rzcl'
        self.parameter['address'] = 'rzcluster.rz.uni-kiel.de:~/'
        self.parameter['queue_system'] = 'SBATCH'
        self.parameter['short_walltime'] = '3:00:00'
        self.parameter['short_batch_class'] = 'express'
        self.parameter['medium_walltime'] = '24:00:00'
        self.parameter['medium_batch_class'] = 'small'
        self.parameter['long_walltime'] = '240:00:00'
        self.parameter['long_batch_class'] = 'long'

    def get_command_line(self):
        """Provides server/cluster command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider a server/cluster.
        """
        new_command_line = str()
        new_command_line += '#' + self.parameter['queue_system'] + ' --output=' \
                            + self.parameter['simulation_directory'] + 'POLARIS.out' + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -j o\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' --time=' \
                            + self.parameter['walltime'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' --mem=' \
                            + self.parameter['ram_usage'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' --nodes=' + self.parameter['node_number'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' --tasks-per-node=1\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' --cpus-per-task=' \
                            + str(self.parameter['nr_threads']) + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' --partition=' + self.parameter[
            'batch_class'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' --job-name=' \
                            + self.parse_args.simulation_name[0:15] + '\n'
        return new_command_line


class AstroServer(Server):
    """This is the server class for the Astro cluster at the ITAP at Kiel university
    """

    def __init__(self, parse_args):
        """Initialisation of the server/cluster parameters.
        """
        Server.__init__(self, parse_args)

        self.parameter['node_name'] = 'ganymed'
        self.parameter['address'] = 'ganymed:~/'
        self.parameter['queue_system'] = 'PBS'
        self.parameter['short_walltime'] = '0:30:00'
        self.parameter['short_batch_class'] = 'wolf-short'
        self.parameter['medium_walltime'] = '3:00:00'
        self.parameter['medium_batch_class'] = 'wolf-medium'
        self.parameter['long_walltime'] = '300:00:00'
        self.parameter['long_batch_class'] = 'wolf-long'
        self.parameter['nr_threads'] = 8

    def get_command_line(self):
        """Provides server/cluster command line for POLARIS .cmd file.

        Returns:
            str: Command line to consider a server/cluster.
        """
        new_command_line = str()
        new_command_line += '#' + self.parameter['queue_system'] + ' -o ' \
                            + self.parameter['simulation_directory'] + 'POLARIS.out' + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -j oe\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -l walltime=' \
                            + self.parameter['walltime'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -l mem=' \
                            + self.parameter['ram_usage'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -l nodes=' + self.parameter['node_number'] \
                            + ':ppn=' + str(self.parameter['nr_threads']) + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -q ' + self.parameter['batch_class'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -N ' \
                            + self.parse_args.simulation_name[0:15] + '\n'
        return new_command_line
