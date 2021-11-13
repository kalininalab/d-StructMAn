class SNV:
    __slots__ = ['new_aa', 'database_id', 'stored', 'tags']

    def __init__(self, new_aa, tags=set(), database_id=None):
        self.new_aa = new_aa
        self.tags = tags.copy()
        self.database_id = database_id
        self.stored = (database_id is not None)

    def fuse(self, other_snv):
        self.tags = self.tags | other_snv.tags
