import matplotlib
if not os.environ.get('DISPLAY'):
    # Use non-interactive Agg backend
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plotCoverage(minorVariantInfo):
    """
    Plot the coverage of a minorVariantInfo instance.
    """
    fig, ax = plt.subplot(1, 1, figsize=(15, 5))

    ax.plot(list(range(len(minorVariantInfo.coveragePerBase))),
            minorVariantInfo.coveragePerBase, '-')

    ax.set_title(minorVariantInfo.name)
    ax.set_xlabel('Position')
    ax.set_ylabel('Coverage')
