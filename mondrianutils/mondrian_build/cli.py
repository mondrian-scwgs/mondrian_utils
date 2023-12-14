import click
import mondrianutils.mondrian_build


@click.group()
def cli():
    pass

@cli.command()
@click.option('--metrics', required=True)
@click.option('--metrics_ref', required=True)
@click.option('--gc_metrics', required=True)
@click.option('--gc_metrics_ref', required=True)
def compare_alignment(metrics, metrics_ref, gc_metrics, gc_metrics_ref):
    mondrianutils.mondrian_build.compare_alignment(
        metrics, metrics_ref,
        gc_metrics, gc_metrics_ref
    )


@cli.command()
@click.option('--reads', required=True)
@click.option('--reads_ref', required=True)
@click.option('--metrics', required=True)
@click.option('--metrics_ref', required=True)
def compare_hmmcopy(reads, reads_ref, metrics, metrics_ref):
    mondrianutils.mondrian_build.compare_hmmcopy(
        reads, reads_ref,
        metrics, metrics_ref
    )


@cli.command()
@click.option('--museq', required=True)
@click.option('--museq_ref', required=True)
@click.option('--mutect', required=True)
@click.option('--mutect_ref', required=True)
@click.option('--strelka_snv', required=True)
@click.option('--strelka_snv_ref', required=True)
@click.option('--strelka_indel', required=True)
@click.option('--strelka_indel_ref', required=True)
def compare_variant_calling(museq, museq_ref, mutect, mutect_ref, strelka_snv, strelka_snv_ref, strelka_indel,
                                strelka_indel_ref):
    mondrianutils.mondrian_build.compare_variant_calling(
        museq, museq_ref,
        mutect, mutect_ref,
        strelka_snv, strelka_snv_ref,
        strelka_indel, strelka_indel_ref
    )


@cli.command()
@click.option('--destruct', required=True)
@click.option('--destruct_ref', required=True)
@click.option('--lumpy', required=True)
@click.option('--lumpy_ref', required=True)
@click.option('--gridss', required=True)
@click.option('--gridss_ref', required=True)
@click.option('--svaba', required=True)
@click.option('--svaba_ref', required=True)
def compare_breakpoint_calling(destruct, destruct_ref, lumpy, lumpy_ref, gridss, gridss_ref, svaba, svaba_ref):
    mondrianutils.mondrian_build.compare_breakpoint_calling(
        destruct, destruct_ref,
        lumpy, lumpy_ref,
        gridss, gridss_ref,
        svaba, svaba_ref,
    )


@cli.command()
@click.option('--genotyper', required=True)
@click.option('--genotyper_ref', required=True)
@click.option('--vartrix', required=True)
@click.option('--vartrix_ref', required=True)
def compare_snv_genotyping(genotyper, genotyper_ref, vartrix, vartrix_ref):
    mondrianutils.mondrian_build.compare_snv_genotyping(
        genotyper, genotyper_ref,
        vartrix, vartrix_ref
    )


@cli.command()
@click.option('--genotyper', required=True)
@click.option('--genotyper_ref', required=True)
def compare_sv_genotyping(genotyper, genotyper_ref):
    mondrianutils.mondrian_build.compare_sv_genotyping(
        genotyper, genotyper_ref
    )


@cli.command()
@click.option('--cells_yaml', required=True)
def compare_normalizer(cells_yaml):
    mondrianutils.mondrian_build.compare_normalizer(
        cells_yaml
    )
