from setuptools import find_packages, setup

setup(
    name="euclid_windows",
    version="1.0",
    description="EUCLID photometric window functions",
    zip_safe=False,
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19",
        "scipy",
	"typing",
        "camb>=1.3.5",
    ],
)
