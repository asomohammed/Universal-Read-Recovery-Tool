from setuptools import setup, find_packages

setup(
    name="universal_read_recovery",
    version="1.0.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    entry_points={
        'console_scripts': [
            'universal-read-recovery=universal_read_recovery.__main__:main',
        ],
    },
    author="Your Name",
    author_email="your.email@example.com",
    description="A tool for recovering unassigned sequencing reads",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/Universal-Read-Recovery-Tool",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)