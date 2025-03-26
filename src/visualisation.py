import matplotlib.pyplot as plt
import os

def plot_results(cropped_map, binary_mask, centroids, output_dir="plots"):
    """Plot the cropped map, binary mask, and centroids."""
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    # Plot the first set of images
    axs[0, 0].imshow(cropped_map.data, cmap='gray')
    axs[0, 0].set_title('Cropped Map')

    axs[0, 1].imshow(binary_mask, cmap='gray')
    axs[0, 1].set_title('Binary Mask')

    axs[0, 2].imshow(binary_mask, cmap='gray')
    axs[0, 2].imshow(cropped_map.data, cmap='gray', alpha=0.5)
    axs[0, 2].set_title('Binary Mask Overlay')

    # Plot the binary mask with centroids
    axs[1, 0].imshow(binary_mask, cmap='gray')
    axs[1, 0].scatter(centroids[:, 1], centroids[:, 0], s=10, c='red', marker='x', linewidths=1)
    axs[1, 0].set_title('Binary Mask with Centroids')

    # Save the plot
    output_path = os.path.join(output_dir, "output_plot.png")
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")