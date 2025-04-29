from setuptools import setup, find_packages

setup(
    name="carbon_easy_builder",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A tool for building carbon-based structures",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/carbon_easy_builder",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
) 