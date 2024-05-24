from setuptools import setup


def read_requirements(path):
    with open(path, "r") as f:
        return [line.strip() for line in f if not line.isspace()]


with open("README.md", "r", encoding="UTF-8") as fh:
    long_description = fh.read()


setup(
    name="romeo",
    version="1.0",
    python_requires=">=3.8",
    install_requires=read_requirements("requirements.txt"),
    packages=["romeo"],
    author="Qianhua Zhu, Hailin Pan",
    author_email="zhuqianhua@genomics.cn",
    description="A Robust Marker Gene Selection by Harmonizing Expression Levels and Positive Ratio in Single-Cell Resolution Transcriptome",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GPL",
    url="https://github.com/BrainStOrmics/Romeo",
)
