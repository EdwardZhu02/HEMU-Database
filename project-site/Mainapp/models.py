# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey and OneToOneField has `on_delete` set to the desired behavior
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from django.db import models


class CoixExp(models.Model):
    gene = models.TextField(blank=True, null=True)
    sample_id = models.TextField(blank=True, null=True)
    fpkm = models.FloatField(blank=True, null=True)
    tpm = models.FloatField(blank=True, null=True)


class CoixSamp(models.Model):
    sample_id = models.TextField(blank=True, null=True)
    sample_download_link = models.TextField(blank=True, null=True)
    sample_tissue = models.TextField(blank=True, null=True)
    sample_project_number = models.TextField(blank=True, null=True)
    sample_description = models.TextField(blank=True, null=True)


class ZeaExp(models.Model):
    gene = models.TextField(blank=True, null=True)
    sample_id = models.TextField(blank=True, null=True)
    fpkm = models.FloatField(blank=True, null=True)
    tpm = models.FloatField(blank=True, null=True)


class ZeaSamp(models.Model):
    sample_id = models.TextField(blank=True, null=True)
    sample_download_link = models.TextField(blank=True, null=True)
    sample_tissue = models.TextField(blank=True, null=True)
    sample_project_number = models.TextField(blank=True, null=True)
    sample_description = models.TextField(blank=True, null=True)
