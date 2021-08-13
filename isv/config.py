import pathlib


class Settings:
    def __init__(self):
        root_dir = str(pathlib.Path(__file__).parent.absolute())
        self.model_dir = root_dir + '/models/'
        self.data_dir = root_dir + '/data/'
        self.valid_chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        self.chromosome_dict = dict(zip(self.valid_chromosomes, range(1, 25)))

        self.gene_type_dict = {'protein_coding': 0, 'pseudogene': 1, 'lncRNA': 2,
                               'miRNA': 3, 'rRNA': 4, 'snRNA': 5, 'unevaluated': 6}

        self.regulatory_type_dict = {'CTCF_binding_site': 7, 'TF_binding_site': 8,
                                     'curated': 9, 'enhancer': 10, 'open_chromatin_region': 11,
                                     'promoter': 12, 'promoter_flanking_region': 13}
        
        self.attributes = [
            'gencode_genes',
            'protein_coding',
            'pseudogenes',
            'mirna',
            'lncrna',
            'rrna',
            'snrna',
            'morbid_genes',
            'disease_associated_genes',
            'hi_genes',  # ONLY FOR LOSSES
            'regions_HI',  # ONLY FOR LOSSES
            'regions_TS',  # ONLY FOR GAINS
            'regulatory',
            'regulatory_enhancer',
            'regulatory_open_chromatin_region',
            'regulatory_promoter',
            'regulatory_promoter_flanking_region',
            'regulatory_ctcf_binding_site',
            'regulatory_tf_binding_site',
            'regulatory_curated'
        ]


settings = Settings()