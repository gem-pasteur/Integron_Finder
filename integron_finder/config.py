import os


from integron_finder import utils

class Config(object):

    def __init__(self, args):
        self._model_len = None  # model_len cache, because it's computation is "heavy" (open file)
        self._args = args

        _prefix_share = '$PREFIXSHARE'
        if 'INTEGRON_HOME' in os.environ and os.environ['INTEGRON_HOME']:
            _prefix_share = os.environ['INTEGRON_HOME']
        elif _prefix_share.endswith('PREFIXSHARE'):
            _prefix_share = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        self._prefix_data = os.path.join(_prefix_share, 'data')

    def __getattr__(self, item):
        try:
            attr = getattr(self._args, item)
            return attr
        except AttributeError:
            raise AttributeError("config object has no attribute '{}'".format(item))

    @property
    def replicon_path(self):
        return os.path.abspath(self._args.replicon)

    @property
    def replicon_name(self):
        in_dir, sequence_file = os.path.split(self.replicon_path)
        replicon_name, extension = os.path.splitext(sequence_file)
        return replicon_name

    @property
    def input_dir(self):
        in_dir, sequence_file = os.path.split(self.replicon_path)
        return in_dir

    @property
    def outdir(self):
        return os.path.abspath(self._args.outdir)

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
        return os.path.join(self._prefix_data, "Models")

    @property
    def model_integrase(self):
        return os.path.join(self.model_dir, "integron_integrase.hmm")

    @property
    def model_phage_int(self):
        return os.path.join(self.model_dir, "phage-int.hmm")

    @property
    def model_attc_path(self):
        try:
            self._args.attc_model
        except AttributeError:
            raise RuntimeError("'model_attc' is not define.")
        if len(self._args.attc_model.split(os.sep)) > 1:  # contain path
            model_attc = self._args.attc_model
        else:
            model_attc = os.path.join(self.model_dir, self._args.attc_model)
        return model_attc

    @property
    def model_attc_name(self):
        try:
            self._args.attc_model
        except AttributeError:
            raise RuntimeError("'model_attc' is not define.")
        return utils.get_name_from_path(self.model_attc_path)

    @property
    def model_len(self):
        try:
            self._args.attc_model
        except AttributeError:
            raise RuntimeError("'model_attc' is not define.")
        if self._model_len is None:
            model_len = utils.model_len(self.model_attc_path)
            self._model_len = model_len
            return model_len
        else:
            return self._model_len

    @property
    def func_annot_path(self):
        return os.path.join(self._prefix_data, "Functional_annotation")
