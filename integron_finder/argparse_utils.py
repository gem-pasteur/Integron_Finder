import integron_finder
import argparse


class MyVersionAction(argparse._VersionAction):

    def __call__(self, parser, namespace, values, option_string=None):
        """
        Override the :meth:`argparse._VersionAction.__call__`
        to allow to call get_version_message with the

        * --hmmsearch
        * --cmsearch
        * --prodigal

        value provided on the command line

        I ensure that the --version option is the last argument
        so all other options are parsed and are in namespace
        then --version is parsed then code below is executed.
        """
        version = integron_finder.get_version_message(hmmsearch=namespace.hmmsearch,
                                                      cmsearch=namespace.cmsearch,
                                                      prodigal=namespace.prodigal)
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), argparse._sys.stdout)
        parser.exit()

