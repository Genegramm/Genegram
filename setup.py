from pathlib import Path

from setuptools import setup, find_packages

root = Path(__file__).parent.resolve()

with open(root / "README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

with open(root / "genegram/__init__.py") as f:
    for line in f:
        if line.startswith("__version__"):
            version = line.strip().split()[-1][1:-1]
            break

name = "genegram"

description = (
    "Project for genomic sequences analysis "
    "employing the combination of formal grammars and neural networks"
)

authors = {
    "vadyushkins": ("Vadim Abzalov", "vadim.i.abzalov@gmail.com"),
    "LuninaPolina": ("Polina Lunina", "lunina_polina@mail.ru"),
}

url = "https://github.com/JetBrains-Research/Genegram"

project_urls = {
    "Documentation": "https://github.com/JetBrains-Research/Genegram",
    "Source Code": "https://github.com/JetBrains-Research/Genegram",
    "Bug Tracker": "https://github.com/JetBrains-Research/Genegram/issues",
}

platforms = ["Linux", "Mac OSX", "Unix"]

keywords = []

classifiers = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Education",
    "License :: OSI Approved :: Apache Software License",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS",
    "Operating System :: Unix",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3 :: Only",
    "Topic :: Software Development :: Libraries :: Python Modules",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Information Analysis",
]


def parse_requirements_file(filename):
    with open(filename) as fid:
        requires = [l.strip() for l in fid.readlines() if not l.startswith("#")]

    return requires


install_requires = parse_requirements_file(root / "requirements" / "default.txt")
extras_require = {
    dep: parse_requirements_file(root / "requirements" / f"{dep}.txt")
    for dep in ["developer", "test"]
}

if __name__ == "__main__":
    setup(
        name=name,
        description=description,
        long_description=long_description,
        author=authors["vadyushkins"][0],
        author_email=authors["vadyushkins"][1],
        maintainer=authors["LuninaPolina"][0],
        maintainer_email=authors["LuninaPolina"][1],
        version=version,
        keywords=keywords,
        packages=find_packages(),
        package_data={"": ["weights/*.h5"]},
        platforms=platforms,
        url=url,
        project_urls=project_urls,
        classifiers=classifiers,
        install_requires=install_requires,
        extras_require=extras_require,
        python_requires=">=3.8",
        zip_safe=False,
    )
