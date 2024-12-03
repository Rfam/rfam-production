# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin sqlcustom [app_label]'
# into your database.
from __future__ import unicode_literals

from django.db import models


class AnnotatedFile(models.Model):
    rfam_acc = models.ForeignKey("Family", db_column="rfam_acc")
    seed = models.TextField()
    cm = models.TextField()
    full = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "_annotated_file"


class FamilyFile(models.Model):
    rfam_acc = models.ForeignKey("Family", db_column="rfam_acc")
    seed = models.TextField()
    cm = models.TextField()

    class Meta:
        managed = False
        db_table = "_family_file"


class Lock(models.Model):
    locked = models.IntegerField()
    locker = models.CharField(max_length=10)
    allowcommits = models.IntegerField(
        db_column="allowCommits"
    )  # Field name made lowercase.
    alsoallow = models.TextField(
        db_column="alsoAllow", blank=True, null=True
    )  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = "_lock"


class Overlap(models.Model):
    auto_overlap = models.AutoField(primary_key=True)
    id = models.CharField(max_length=40, blank=True, null=True)
    description = models.CharField(max_length=100, blank=True, null=True)
    author = models.TextField(blank=True, null=True)
    comment = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "_overlap"


class OverlapMembership(models.Model):
    rfam_acc = models.ForeignKey("Family", db_column="rfam_acc")
    auto_overlap = models.ForeignKey(Overlap, db_column="auto_overlap")

    class Meta:
        managed = False
        db_table = "_overlap_membership"


class PostProcess(models.Model):
    rfam_acc = models.ForeignKey("Family", db_column="rfam_acc")
    author = models.CharField(max_length=45)
    uuid = models.CharField(max_length=45)
    status = models.CharField(max_length=4)
    created = models.DateTimeField()
    opened = models.DateTimeField(blank=True, null=True)
    closed = models.DateTimeField(blank=True, null=True)
    message = models.TextField(blank=True, null=True)
    lsf_id = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "_post_process"


class AlignmentAndTree(models.Model):
    rfam_acc = models.ForeignKey("Family", db_column="rfam_acc")
    type = models.CharField(max_length=9)
    alignment = models.TextField(blank=True, null=True)
    tree = models.TextField(blank=True, null=True)
    treemethod = models.TextField(blank=True, null=True)
    average_length = models.FloatField(blank=True, null=True)
    percent_id = models.FloatField(blank=True, null=True)
    number_of_sequences = models.BigIntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "alignment_and_tree"


class AuthGroup(models.Model):
    name = models.CharField(unique=True, max_length=80)

    class Meta:
        managed = False
        db_table = "auth_group"


class AuthGroupPermissions(models.Model):
    group_id = models.ForeignKey(AuthGroup)
    permission_id = models.ForeignKey("AuthPermission")

    class Meta:
        managed = False
        db_table = "auth_group_permissions"
        unique_together = (("group_id", "permission_id"),)


class AuthPermission(models.Model):
    name = models.CharField(max_length=255)
    content_type_id = models.ForeignKey("DjangoContentType")
    codename = models.CharField(max_length=100)

    class Meta:
        managed = False
        db_table = "auth_permission"
        unique_together = (("content_type_id", "codename"),)


class AuthUser(models.Model):
    password = models.CharField(max_length=128)
    last_login = models.DateTimeField(blank=True, null=True)
    is_superuser = models.IntegerField()
    username = models.CharField(unique=True, max_length=30)
    first_name = models.CharField(max_length=30)
    last_name = models.CharField(max_length=30)
    email = models.CharField(max_length=254)
    is_staff = models.IntegerField()
    is_active = models.IntegerField()
    date_joined = models.DateTimeField()

    class Meta:
        managed = False
        db_table = "auth_user"


class AuthUserGroups(models.Model):
    user_id = models.ForeignKey(AuthUser)
    group_id = models.ForeignKey(AuthGroup)

    class Meta:
        managed = False
        db_table = "auth_user_groups"
        unique_together = (("user_id", "group_id"),)


class AuthUserUserPermissions(models.Model):
    user_id = models.ForeignKey(AuthUser)
    permission_id = models.ForeignKey(AuthPermission)

    class Meta:
        managed = False
        db_table = "auth_user_user_permissions"
        unique_together = (("user_id", "permission_id"),)


class Clan(models.Model):
    clan_acc = models.CharField(primary_key=True, max_length=7)
    id = models.CharField(max_length=40, blank=True, null=True)
    previous_id = models.TextField(blank=True, null=True)
    description = models.CharField(max_length=100, blank=True, null=True)
    author = models.TextField(blank=True, null=True)
    comment = models.TextField(blank=True, null=True)
    created = models.DateTimeField()
    updated = models.DateTimeField()

    class Meta:
        managed = False
        db_table = "clan"


class ClanDatabaseLink(models.Model):
    clan_acc = models.ForeignKey(Clan, db_column="clan_acc")
    db_id = models.TextField()
    comment = models.TextField(blank=True, null=True)
    db_link = models.TextField()
    other_params = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "clan_database_link"


class ClanLiteratureReference(models.Model):
    clan_acc = models.ForeignKey(Clan, db_column="clan_acc")
    pmid = models.ForeignKey("LiteratureReference", db_column="pmid")
    comment = models.TextField(blank=True, null=True)
    order_added = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "clan_literature_reference"


class ClanMembership(models.Model):
    clan_acc = models.ForeignKey(Clan, db_column="clan_acc")
    rfam_acc = models.ForeignKey("Family", db_column="rfam_acc", unique=True)

    class Meta:
        managed = False
        db_table = "clan_membership"


class DatabaseLink(models.Model):
    rfam_acc = models.ForeignKey("Family", db_column="rfam_acc")
    db_id = models.TextField()
    comment = models.TextField(blank=True, null=True)
    db_link = models.TextField()
    other_params = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "database_link"


class DbVersion(models.Model):
    rfam_release = models.FloatField(primary_key=True)
    rfam_release_date = models.DateTimeField()
    number_families = models.IntegerField()
    embl_release = models.TextField()
    genome_collection_date = models.DateTimeField(blank=True, null=True)
    refseq_version = models.IntegerField(blank=True, null=True)
    pdb_date = models.DateTimeField(blank=True, null=True)
    infernal_version = models.CharField(max_length=45, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "db_version"


class DeadClan(models.Model):
    clan_acc = models.CharField(unique=True, max_length=7)
    clan_id = models.CharField(max_length=40)
    comment = models.TextField(blank=True, null=True)
    forward_to = models.CharField(max_length=7, blank=True, null=True)
    user = models.TextField()

    class Meta:
        managed = False
        db_table = "dead_clan"


class DeadFamily(models.Model):
    rfam_acc = models.CharField(unique=True, max_length=7)
    rfam_id = models.CharField(max_length=40)
    comment = models.TextField(blank=True, null=True)
    forward_to = models.CharField(max_length=7, blank=True, null=True)
    title = models.CharField(max_length=150, blank=True, null=True)
    user = models.TextField()

    class Meta:
        managed = False
        db_table = "dead_family"


class DjangoAdminLog(models.Model):
    action_time = models.DateTimeField()
    object_id = models.TextField(blank=True, null=True)
    object_repr = models.CharField(max_length=200)
    action_flag = models.SmallIntegerField()
    change_message = models.TextField()
    content_type = models.ForeignKey("DjangoContentType", blank=True, null=True)
    user = models.ForeignKey(AuthUser)

    class Meta:
        managed = False
        db_table = "django_admin_log"


class DjangoContentType(models.Model):
    app_label = models.CharField(max_length=100)
    model = models.CharField(max_length=100)

    class Meta:
        managed = False
        db_table = "django_content_type"
        unique_together = (("app_label", "model"),)


class DjangoMigrations(models.Model):
    app = models.CharField(max_length=255)
    name = models.CharField(max_length=255)
    applied = models.DateTimeField()

    class Meta:
        managed = False
        db_table = "django_migrations"


class DjangoSession(models.Model):
    session_key = models.CharField(primary_key=True, max_length=40)
    session_data = models.TextField()
    expire_date = models.DateTimeField()

    class Meta:
        managed = False
        db_table = "django_session"


class Family(models.Model):
    rfam_acc = models.CharField(primary_key=True, max_length=7)
    rfam_id = models.CharField(max_length=40)
    auto_wiki = models.ForeignKey("Wikitext", db_column="auto_wiki")
    description = models.CharField(max_length=75, blank=True, null=True)
    author = models.TextField(blank=True, null=True)
    seed_source = models.TextField(blank=True, null=True)
    gathering_cutoff = models.FloatField(blank=True, null=True)
    trusted_cutoff = models.FloatField(blank=True, null=True)
    noise_cutoff = models.FloatField(blank=True, null=True)
    comment = models.TextField(blank=True, null=True)
    previous_id = models.TextField(blank=True, null=True)
    cmbuild = models.TextField(blank=True, null=True)
    cmcalibrate = models.TextField(blank=True, null=True)
    cmsearch = models.TextField(blank=True, null=True)
    num_seed = models.BigIntegerField(blank=True, null=True)
    num_full = models.BigIntegerField(blank=True, null=True)
    num_genome_seq = models.BigIntegerField(blank=True, null=True)
    num_refseq = models.BigIntegerField(blank=True, null=True)
    type = models.CharField(max_length=50, blank=True, null=True)
    structure_source = models.TextField(blank=True, null=True)
    number_of_species = models.BigIntegerField(blank=True, null=True)
    number_3d_structures = models.IntegerField(blank=True, null=True)
    tax_seed = models.TextField(blank=True, null=True)
    ecmli_lambda = models.FloatField(blank=True, null=True)
    ecmli_mu = models.FloatField(blank=True, null=True)
    ecmli_cal_db = models.IntegerField(blank=True, null=True)
    ecmli_cal_hits = models.IntegerField(blank=True, null=True)
    maxl = models.IntegerField(blank=True, null=True)
    clen = models.IntegerField(blank=True, null=True)
    match_pair_node = models.IntegerField(blank=True, null=True)
    hmm_tau = models.FloatField(blank=True, null=True)
    hmm_lambda = models.FloatField(blank=True, null=True)
    created = models.DateTimeField()
    updated = models.DateTimeField()

    class Meta:
        managed = False
        db_table = "family"


class FamilyLiteratureReference(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    pmid = models.ForeignKey("LiteratureReference", db_column="pmid")
    comment = models.TextField(blank=True, null=True)
    order_added = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "family_literature_reference"


class FamilyLong(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    referenece_structure = models.TextField(blank=True, null=True)
    reference_sequence = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "family_long"


class FamilyNcbi(models.Model):
    ncbi = models.ForeignKey("Taxonomy")
    rfam_id = models.CharField(max_length=40, blank=True, null=True)
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")

    class Meta:
        managed = False
        db_table = "family_ncbi"


class Features(models.Model):
    rfamseq_acc = models.ForeignKey("Rfamseq", db_column="rfamseq_acc")
    database_id = models.CharField(max_length=50)
    primary_id = models.CharField(max_length=100)
    secondary_id = models.CharField(max_length=255, blank=True, null=True)
    feat_orient = models.IntegerField()
    feat_start = models.BigIntegerField()
    feat_end = models.BigIntegerField()
    quaternary_id = models.CharField(max_length=150, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "features"


class FullRegion(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    rfamseq_acc = models.ForeignKey("Rfamseq", db_column="rfamseq_acc")
    seq_start = models.BigIntegerField()
    seq_end = models.BigIntegerField()
    bit_score = models.FloatField()
    evalue_score = models.CharField(max_length=15)
    cm_start = models.IntegerField()
    cm_end = models.IntegerField()
    truncated = models.CharField(max_length=2)
    type = models.CharField(max_length=4)
    is_significant = models.IntegerField()

    class Meta:
        managed = False
        db_table = "full_region"


class Genome(models.Model):
    upid = models.CharField(primary_key=True, max_length=20)
    assembly_acc = models.CharField(max_length=20)
    assembly_version = models.IntegerField(blank=True, null=True)
    wgs_acc = models.CharField(max_length=20, blank=True, null=True)
    wgs_version = models.IntegerField(blank=True, null=True)
    assembly_name = models.CharField(max_length=100, blank=True, null=True)
    assembly_level = models.CharField(max_length=15, blank=True, null=True)
    study_ref = models.CharField(max_length=20, blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    total_length = models.BigIntegerField(blank=True, null=True)
    ungapped_length = models.BigIntegerField(blank=True, null=True)
    circular = models.IntegerField(blank=True, null=True)
    ncbi = models.ForeignKey("Taxonomy")
    scientific_name = models.CharField(max_length=100, blank=True, null=True)
    common_name = models.CharField(max_length=200, blank=True, null=True)
    kingdom = models.CharField(max_length=50, blank=True, null=True)
    num_rfam_regions = models.IntegerField(blank=True, null=True)
    num_families = models.IntegerField(blank=True, null=True)
    created = models.DateTimeField()
    updated = models.DateTimeField()

    class Meta:
        managed = False
        db_table = "genome"


class Genseq(models.Model):
    rfamseq_acc = models.CharField(primary_key=True, max_length=20)
    upid = models.CharField(max_length=20)
    chromosome_name = models.CharField(max_length=20)
    chromosome_type = models.CharField(max_length=20)

    class Meta:
        managed = False
        db_table = "genseq"


class HtmlAlignment(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    type = models.CharField(max_length=16)
    html = models.TextField(blank=True, null=True)
    block = models.IntegerField()
    html_alignmentscol = models.CharField(max_length=45, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "html_alignment"


class Keywords(models.Model):
    rfam_acc = models.CharField(primary_key=True, max_length=7)
    rfam_id = models.CharField(max_length=40, blank=True, null=True)
    description = models.CharField(max_length=100, blank=True, null=True)
    rfam_general = models.TextField(blank=True, null=True)
    literature = models.TextField(blank=True, null=True)
    wiki = models.TextField(blank=True, null=True)
    pdb_mappings = models.TextField(blank=True, null=True)
    clan_info = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "keywords"


class LiteratureReference(models.Model):
    pmid = models.AutoField(primary_key=True)
    title = models.TextField(blank=True, null=True)
    author = models.TextField(blank=True, null=True)
    journal = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "literature_reference"


class MatchesAndFasta(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    match_list = models.TextField(blank=True, null=True)
    fasta = models.TextField(blank=True, null=True)
    type = models.CharField(max_length=7)

    class Meta:
        managed = False
        db_table = "matches_and_fasta"


class Motif(models.Model):
    motif_acc = models.CharField(primary_key=True, max_length=7)
    motif_id = models.CharField(max_length=40, blank=True, null=True)
    description = models.CharField(max_length=75, blank=True, null=True)
    author = models.TextField(blank=True, null=True)
    seed_source = models.TextField(blank=True, null=True)
    gathering_cutoff = models.FloatField(blank=True, null=True)
    trusted_cutoff = models.FloatField(blank=True, null=True)
    noise_cutoff = models.FloatField(blank=True, null=True)
    cmbuild = models.TextField(blank=True, null=True)
    cmcalibrate = models.TextField(blank=True, null=True)
    type = models.CharField(max_length=50, blank=True, null=True)
    num_seed = models.BigIntegerField(blank=True, null=True)
    average_id = models.FloatField(blank=True, null=True)
    average_sqlen = models.FloatField(blank=True, null=True)
    ecmli_lambda = models.FloatField(blank=True, null=True)
    ecmli_mu = models.FloatField(blank=True, null=True)
    ecmli_cal_db = models.IntegerField(blank=True, null=True)
    ecmli_cal_hits = models.IntegerField(blank=True, null=True)
    maxl = models.IntegerField(blank=True, null=True)
    clen = models.IntegerField(blank=True, null=True)
    match_pair_node = models.IntegerField(blank=True, null=True)
    hmm_tau = models.FloatField(blank=True, null=True)
    hmm_lambda = models.FloatField(blank=True, null=True)
    wiki = models.CharField(max_length=80, blank=True, null=True)
    created = models.DateTimeField()
    updated = models.DateTimeField()

    class Meta:
        managed = False
        db_table = "motif"


class MotifDatabaseLink(models.Model):
    motif_acc = models.ForeignKey(Motif, db_column="motif_acc")
    db_id = models.TextField()
    comment = models.TextField(blank=True, null=True)
    db_link = models.TextField()
    other_params = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "motif_database_link"


class MotifFamilyStats(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    motif_acc = models.ForeignKey("MotifOld", db_column="motif_acc")
    num_hits = models.IntegerField(blank=True, null=True)
    frac_hits = models.DecimalField(
        max_digits=4, decimal_places=3, blank=True, null=True
    )
    sum_bits = models.DecimalField(
        max_digits=12, decimal_places=3, blank=True, null=True
    )
    avg_weight_bits = models.DecimalField(
        max_digits=12, decimal_places=3, blank=True, null=True
    )

    class Meta:
        managed = False
        db_table = "motif_family_stats"


class MotifFile(models.Model):
    motif_acc = models.ForeignKey(Motif, db_column="motif_acc")
    seed = models.TextField()
    cm = models.TextField()

    class Meta:
        managed = False
        db_table = "motif_file"


class MotifLiterature(models.Model):
    motif_acc = models.ForeignKey("MotifOld", db_column="motif_acc")
    pmid = models.ForeignKey(LiteratureReference, db_column="pmid")
    comment = models.TextField(blank=True, null=True)
    order_added = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "motif_literature"


class MotifMatches(models.Model):
    motif_acc = models.ForeignKey("MotifOld", db_column="motif_acc")
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    rfamseq_acc = models.ForeignKey("Rfamseq", db_column="rfamseq_acc")
    rfamseq_start = models.BigIntegerField(blank=True, null=True)
    rfamseq_stop = models.BigIntegerField(blank=True, null=True)
    query_start = models.IntegerField(blank=True, null=True)
    query_stop = models.IntegerField(blank=True, null=True)
    motif_start = models.IntegerField(blank=True, null=True)
    motif_stop = models.IntegerField(blank=True, null=True)
    e_value = models.CharField(max_length=15, blank=True, null=True)
    bit_score = models.FloatField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "motif_matches"


class MotifOld(models.Model):
    motif_acc = models.CharField(primary_key=True, max_length=7)
    motif_id = models.CharField(max_length=40, blank=True, null=True)
    description = models.CharField(max_length=75, blank=True, null=True)
    author = models.TextField(blank=True, null=True)
    seed_source = models.TextField(blank=True, null=True)
    gathering_cutoff = models.FloatField(blank=True, null=True)
    trusted_cutoff = models.FloatField(blank=True, null=True)
    noise_cutoff = models.FloatField(blank=True, null=True)
    cmbuild = models.TextField(blank=True, null=True)
    cmcalibrate = models.TextField(blank=True, null=True)
    type = models.CharField(max_length=50, blank=True, null=True)
    ecmli_lambda = models.FloatField(blank=True, null=True)
    ecmli_mu = models.FloatField(blank=True, null=True)
    ecmli_cal_db = models.IntegerField(blank=True, null=True)
    ecmli_cal_hits = models.IntegerField(blank=True, null=True)
    maxl = models.IntegerField(blank=True, null=True)
    clen = models.IntegerField(blank=True, null=True)
    match_pair_node = models.IntegerField(blank=True, null=True)
    hmm_tau = models.FloatField(blank=True, null=True)
    hmm_lambda = models.FloatField(blank=True, null=True)
    created = models.DateTimeField()
    updated = models.DateTimeField()

    class Meta:
        managed = False
        db_table = "motif_old"


class MotifPdb(models.Model):
    motif_acc = models.ForeignKey(MotifOld, db_column="motif_acc")
    pdb_id = models.CharField(max_length=4)
    chain = models.CharField(max_length=4, blank=True, null=True)
    pdb_start = models.IntegerField(blank=True, null=True)
    pdb_end = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "motif_pdb"


class MotifSsImage(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    motif_acc = models.ForeignKey(MotifOld, db_column="motif_acc")
    image = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "motif_ss_image"


class Pdb(models.Model):
    pdb_id = models.CharField(primary_key=True, max_length=4)
    keywords = models.TextField(blank=True, null=True)
    title = models.TextField(blank=True, null=True)
    date = models.TextField(blank=True, null=True)
    resolution = models.DecimalField(
        max_digits=5, decimal_places=2, blank=True, null=True
    )
    method = models.TextField(blank=True, null=True)
    author = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "pdb"


class PdbFullRegion(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    pdb_id = models.CharField(max_length=4)
    chain = models.CharField(max_length=4, blank=True, null=True)
    pdb_start = models.IntegerField()
    pdb_end = models.IntegerField()
    bit_score = models.FloatField()
    evalue_score = models.CharField(max_length=15)
    cm_start = models.IntegerField()
    cm_end = models.IntegerField()
    hex_colour = models.CharField(max_length=6, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "pdb_full_region"


class PdbRfamReg(models.Model):
    auto_pdb_reg = models.AutoField(primary_key=True)
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    pdb_seq = models.ForeignKey("PdbSequence", db_column="pdb_seq")
    pdb = models.ForeignKey(Pdb)
    chain = models.CharField(max_length=4, blank=True, null=True)
    pdb_res_start = models.IntegerField(blank=True, null=True)
    pdb_res_end = models.IntegerField(blank=True, null=True)
    rfamseq_acc = models.ForeignKey("Rfamseq", db_column="rfamseq_acc")
    seq_start = models.BigIntegerField(blank=True, null=True)
    seq_end = models.BigIntegerField(blank=True, null=True)
    hex_colour = models.CharField(max_length=6, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "pdb_rfam_reg"


class PdbSequence(models.Model):
    pdb_seq = models.CharField(primary_key=True, max_length=6)
    pdb = models.ForeignKey(Pdb)
    chain = models.CharField(max_length=1, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "pdb_sequence"


class ProcessedData(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    cm = models.TextField(blank=True, null=True)
    ss_stats_pbp = models.TextField(blank=True, null=True)
    ss_stats_seq = models.TextField(blank=True, null=True)
    ss_stats_fam = models.TextField(blank=True, null=True)
    scores_graph = models.TextField(blank=True, null=True)
    genome_full = models.TextField(blank=True, null=True)
    genome_full_md5 = models.CharField(max_length=32, blank=True, null=True)
    refseq_full = models.TextField(blank=True, null=True)
    refseq_full_md5 = models.CharField(max_length=32, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "processed_data"


class Refseq(models.Model):
    refseq_acc = models.CharField(primary_key=True, max_length=14)
    description = models.TextField(blank=True, null=True)
    species = models.TextField(blank=True, null=True)
    ncbi_taxid = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "refseq"


class RefseqFullRegion(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    refseq_acc = models.ForeignKey(Refseq, db_column="refseq_acc")
    seq_start = models.BigIntegerField()
    seq_end = models.BigIntegerField()
    bit_score = models.FloatField()
    evalue_score = models.CharField(max_length=15)
    cm_start = models.IntegerField()
    cm_end = models.IntegerField()
    truncated = models.CharField(max_length=2)

    class Meta:
        managed = False
        db_table = "refseq_full_region"


class Rfamseq(models.Model):
    rfamseq_acc = models.CharField(primary_key=True, max_length=20)
    accession = models.CharField(max_length=15)
    version = models.IntegerField()
    ncbi = models.ForeignKey("Taxonomy")
    mol_type = models.CharField(max_length=15)
    length = models.IntegerField(blank=True, null=True)
    description = models.CharField(max_length=250)
    previous_acc = models.TextField(blank=True, null=True)
    source = models.CharField(max_length=20)

    class Meta:
        managed = False
        db_table = "rfamseq"


class SecondaryStructureImage(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    type = models.CharField(max_length=8, blank=True, null=True)
    image = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = "secondary_structure_image"


class SeedRegion(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    rfamseq_acc = models.ForeignKey(Rfamseq, db_column="rfamseq_acc")
    seq_start = models.BigIntegerField()
    seq_end = models.BigIntegerField()

    class Meta:
        managed = False
        db_table = "seed_region"


class Sunburst(models.Model):
    rfam_acc = models.ForeignKey(Family, db_column="rfam_acc")
    data = models.TextField()
    type = models.CharField(max_length=7)

    class Meta:
        managed = False
        db_table = "sunburst"


class TaxonomicTree(models.Model):
    ncbi_code = models.IntegerField(primary_key=True)
    species = models.CharField(max_length=100, blank=True, null=True)
    taxonomy = models.TextField(blank=True, null=True)
    lft = models.IntegerField(blank=True, null=True)
    rgt = models.IntegerField(blank=True, null=True)
    parent = models.CharField(max_length=100, blank=True, null=True)
    level = models.CharField(max_length=100, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "taxonomic_tree"


class Taxonomy(models.Model):
    ncbi_id = models.IntegerField(primary_key=True)
    species = models.CharField(max_length=100)
    tax_string = models.TextField(blank=True, null=True)
    tree_display_name = models.CharField(max_length=100, blank=True, null=True)
    align_display_name = models.CharField(max_length=50, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "taxonomy"


class TaxonomyWebsearch(models.Model):
    ncbi_id = models.IntegerField(blank=True, null=True)
    species = models.CharField(max_length=100, blank=True, null=True)
    taxonomy = models.TextField(blank=True, null=True)
    lft = models.IntegerField(blank=True, null=True)
    rgt = models.IntegerField(blank=True, null=True)
    parent = models.IntegerField(blank=True, null=True)
    level = models.CharField(max_length=200, blank=True, null=True)
    minimal = models.IntegerField()
    rank = models.CharField(max_length=100, blank=True, null=True)

    class Meta:
        managed = False
        db_table = "taxonomy_websearch"


class Wikitext(models.Model):
    auto_wiki = models.AutoField(primary_key=True)
    title = models.CharField(unique=True, max_length=150)

    class Meta:
        managed = False
        db_table = "wikitext"
