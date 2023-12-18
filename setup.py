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
            'alignment_utils = mondrianutils.alignment.cli:cli',
            'hmmcopy_utils = mondrianutils.hmmcopy.cli:cli',
            'breakpoint_utils = mondrianutils.breakpoint_calling.cli:cli',
            'haplotype_utils = mondrianutils.haplotypes.cli:cli',
            'variant_utils = mondrianutils.variant_calling.cli:cli',
            'snv_genotyping_utils = mondrianutils.snv_genotyping.cli:cli',
            'sv_genotyping_utils = mondrianutils.sv_genotyping.cli:cli',
            'normalizer_utils = mondrianutils.normalizer.cli:cli',
            'io_utils = mondrianutils.io.cli:cli',
            'reference_utils = mondrianutils.reference.cli:cli',
            'dlp_utils = mondrianutils.dlp_utils.cli:cli',
            'mondrian_build_utils = mondrianutils.mondrian_build.cli:cli'
        ]
    },
    package_data={'': ['*.py', '*.R', '*.npz', "*.yaml", "data/*", "*.sh"]}
)
