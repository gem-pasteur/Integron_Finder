import os




class Config(object):

    def __init__(self, args):
        self._args = args

        _prefix_share = '$PREFIXSHARE'
        if 'INTEGRON_HOME' in os.environ and os.environ['INTEGRON_HOME']:
            _prefix_share = os.environ['INTEGRON_HOME']

        self._prefix_data = os.path.join(_prefix_share, 'data')

    def __getattr__(self, item):
        return getattr(self.args, item)

    @property
    def replicon_path(self):
        os.path.abspath(self._args.replicon)

    @property
    def replicon_name(self):
        in_dir, sequence_file = os.path.split(self.replicon_path)
        replicon_name, extension = os.path.splitext(sequence_file)
        return replicon_name

    @property
    def input_dir(self):
        in_dir, sequence_file = os.path.abspath(os.path.split(self.replicon_path))
        return in_dir

    @property
    def mode_name(self):
        if self._args.eagle_eyes or self._args.local_max:
            return "local_max"
        else:
            return "default"

    @property
    def default_topology(self):
        if self._args.circular:
            return 'circ'
        elif self._args.linear:
            return 'lin'
        else:
            return None

    @property
    def model_dir(self):
        return os.path.join(self._prefix_data, "Models/")

    @property
    def model_integrase(self):
        return os.path.join(self.model_dir, "integron_integrase.hmm")

    @property
    def model_phage_int(self):
        return os.path.join(self.model_dir, "phage-int.hmm")

    @property
    def model_attc(self):
        if len(self._args.attc_model.split("/")) > 1: #contain path
            model_attc = self._args.attc_model
        else:
            model_attc = os.path.join(self.model_dir, self._args.attc_model)
        return model_attc

    @property
    def func_annot_path(self):
        os.path.join(self._prefix_data, "Functional_annotation")
