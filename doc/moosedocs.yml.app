###############################################################
# Configuration for creating STORK application documentation. #
###############################################################

# Executable directory or executable
app: ../stork-opt

defaults:

    # Install markdown directory
    install: content

    # The default repository for linking
    details: details

# Directories/applications to include
include:
    stork:
        source: ..
        links:
            Tests:
                - ../test/tests
            Examples:
                - ../examples
