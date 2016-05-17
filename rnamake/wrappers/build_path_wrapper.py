from rnamake import settings
import wrapper

class BuildPathWrapper(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/path_builder"
        super(self.__class__, self).__init__(program_path)

        self.add_cmd_option("mg", "", required=True)


class BuildPathWrapperNew(wrapper.Wrapper):
    def __init__(self):
        program_path = settings.LIB_PATH + "/lib/RNAMake/cmake/build/path_builder_new"
        super(self.__class__, self).__init__(program_path)

        self.add_cmd_option("mg", "", required=True)
        self.add_cmd_option("write_pdbs", False, required=False)
