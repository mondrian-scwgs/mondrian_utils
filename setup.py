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
            'variant_utils = mondrianutils.variant_calling.utils:utils',
            'breakpoint_utils = mondrianutils.breakpoint_calling.utils:utils',
            'alignment_utils = mondrianutils.alignment.utils:utils',
            'hmmcopy_utils = mondrianutils.hmmcopy.utils:utils',
            'csverve_utils = mondrianutils.io.csverve:utils',
            'pdf_utils = mondrianutils.io.pdf:utils',
            'vcf_utils = mondrianutils.io.vcf:utils',
            'bam_utils = mondrianutils.io.bam:utils',
            'haplotype_utils = mondrianutils.haplotypes.utils:utils',
            'snv_genotyping_utils = mondrianutils.snv_genotyping.utils:utils',
            'reference_utils = mondrianutils.reference.utils:utils',
            'dlp_utils = mondrianutils.dlp_utils.utils:utils',
        ]
    },
    package_data={'': ['*.py', '*.R', '*.npz', "*.yaml", "data/*", "*.sh"]}
)
