from pydantic import BaseModel


class FeaturesProps(BaseModel):
    max_grid_height: int = 80
    rec_items_height: int = 12
    plus_strand_level: int = 0
    minus_strand_level: int = 55
    amplicon_pool1_level: int = 0
    amplicon_pool2_level: int = 15
    amplicon_offset: int = 80
    gene_feature_padding: int = 3
