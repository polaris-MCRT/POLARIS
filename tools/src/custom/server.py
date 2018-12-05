#!/usr/bin/env python
# -*- coding: utf-8 -*-

from modules.base import Server

"""Add your defined classes to this dictionary with a unique name
 to use it with PolarisTools.
"""


def update_server_dict(dictionary):
    server_dict = {
        'herschel': HerschelServer,
        'calculus': CalculusServer,
        'cea': CeaServer,
        'custom': CustomServer,
    }
    dictionary.update(server_dict)


class CustomServer(Server):
    """Change this to the server you want to use.
    """

    def __init__(self, parse_args):
        """Initialisation of the server/cluster parameters.
        """
        Server.__init__(self, parse_args)

        # To obtain the 'node_name', enter the following commands on the server/cluster:
        # python
        # import platform
        # platform.node()
        # -> 'custom_system' (e.g. 'nesh-fe2' -> use 'nesh' as 'node_name')
        self.parameter['node_name'] = 'custom_system'
        # This directory is defined as a relative path to your home dir on the server
        self.parameter['server_polaris_dir'] = 'the_polaris_directory'
        self.parameter['address'] = 'custom_server.de:~/'
        self.parameter['queue_system'] = 'PBS'
        self.parameter['short_walltime'] = '2:00:00'
        self.parameter['short_batch_class'] = 'clexpress'
        self.parameter['medium_walltime'] = '48:00:00'
        self.parameter['medium_batch_class'] = 'clmedium'
        self.parameter['long_walltime'] = '100:00:00'
        self.parameter['long_batch_class'] = 'cllong'
        self.parameter['nr_threads'] = 16

    def get_command_line(self):
        """Provides server/cluster command line for POLARIS .cmd file. The syntax depends on the used queue system.

        Returns:
            str: Command line to consider a server/cluster.
        """
        new_command_line = str()
        new_command_line += '#' + self.parameter['queue_system'] + ' -o ' \
                            + self.parameter['simulation_directory'] + \
            'POLARIS.out' + '\n'
        new_command_line += '#' + self.parameter['queue_system'] + ' -j o\n'
        new_command_line += '#' + self.parameter['queue_system'] \
                            + ' -l elapstim_req=' + \
            self.parameter['walltime'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] \
                            + ' -l memsz_job=' + \
            self.parameter['ram_usage'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] \
                            + ' -b ' + self.parameter['node_number'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] \
                            + ' -l cpunum_job=' + \
            str(self.parameter['nr_threads']) + '\n'
        new_command_line += '#' + self.parameter['queue_system'] \
                            + ' -q ' + self.parameter['batch_class'] + '\n'
        new_command_line += '#' + self.parameter['queue_system'] \
                            + ' -N ' + \
            self.parse_args.simulation_name[0:15] + '\n'
        return new_command_line


class CeaServer(Server):
    """This is the server class for the Astro cluster at the ITAP at Kiel university
    """

    def __init__(self, parse_args):
        """Initialisation of the server/cluster parameters.
        """
        Server.__init__(self, parse_args)

        self.parameter['node_name'] = 'sappcw5'
        self.parameter['address'] = 'cea_server:~/'
        self.parameter['server_polaris_dir'] = 'astrophysics/polaris/'


class HerschelServer(Server):
    """This is the server class for the herschel server at IAS Orsay
    """

    def __init__(self, parse_args):
        """Initialisation of the server/cluster parameters.
        """
        Server.__init__(self, parse_args)

        self.parameter['node_name'] = 'glx-herschel'
        self.parameter['address'] = 'herschel_proxy:~/'
        self.parameter['server_polaris_dir'] = 'polaris/'
        self.parameter['queue_system'] = None


class CalculusServer(Server):
    """This is the server class for the calculus server at IAS Orsay
    """

    def __init__(self, parse_args):
        """Initialisation of the server/cluster parameters.
        """
        Server.__init__(self, parse_args)

        self.parameter['node_name'] = 'glx-calcul1'
        self.parameter['address'] = 'glx-calcul1.ias.u-psud.fr:~/'
        self.parameter['server_polaris_dir'] = 'polaris/'
        self.parameter['queue_system'] = None
