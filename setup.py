import versioneer
from setuptools import setup, find_packages

setup(
    name='mondrianutils',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='python utilities for mondrian',
    author='Diljot Grewal',
    author_email='diljot.grewal@gmail.com',
    entry_points={
        'console_scripts': [
            'mondrianutils = mondrianutils.run:main',
            'variant_utils = mondrianutils.variant_calling.utils:cli',
            'breakpoint_utils = mondrianutils.breakpoint_calling.utils:cli',
            'alignment_utils = mondrianutils.alignment.utils:cli',
            'hmmcopy_utils = mondrianutils.hmmcopy.utils:cli',
            'csverve_utils = mondrianutils.io.csverve:cli',
            'pdf_utils = mondrianutils.io.pdf:cli',
            'vcf_utils = mondrianutils.io.vcf:cli',
            'bam_utils = mondrianutils.io.bam:cli',
            'haplotype_utils = mondrianutils.haplotypes.utils:cli',
            'snv_genotyping_utils = mondrianutils.snv_genotyping.utils:cli',
            'sv_genotyping_utils = mondrianutils.sv_genotyping.utils:cli',
            'reference_utils = mondrianutils.reference.utils:cli',
            'dlp_utils = mondrianutils.dlp_utils.utils:cli',
            'normalizer_utils = mondrianutils.normalizer.utils:cli',
            'mondrian_build_utils = mondrianutils.mondrian_build.utils:cli'
        ]
    },
    package_data={'': ['*.py', '*.R', '*.npz', "*.yaml", "data/*", "*.sh"]}
)
