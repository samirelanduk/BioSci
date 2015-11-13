from distutils.core import setup

setup(name="biosci",
      version="0.1.0",
      description="Biological utilities",
      long_description="Assorted tools for dealing with biological information and tools.",
      url="https://github.com/samirelanduk/BioSci",
      author="Sam Ireland",
      author_email="sam.ireland.uk@gmail.com",
      classifiers=["Development Status :: 4 - Beta",
                   "Programming Language :: Python :: 2.6",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3"],
      packages=["biosci"],
      install_requires=["requests"]) #Soon anyway
