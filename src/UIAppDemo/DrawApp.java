package UIAppDemo;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.event.*;
import MLComp.MLOps;
import MLComp.Neuron;
import MLComp.ActivationFunctions_Folder.ActivationFunctions;

public class DrawApp extends JFrame {

    private DrawingPanel drawingPanel;
    private JTextArea resultArea;
    private Neuron[][] network;
    private static final String MODEL_PATH = "Models//Optimus.ser";

    public DrawApp() {
        setTitle("Neural Network Digit Drawer");
        setSize(600, 400);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());

        // Initialize components
        drawingPanel = new DrawingPanel();
        JPanel controlPanel = new JPanel();
        JButton clearButton = new JButton("Clear");
        JButton predictButton = new JButton("Predict");
        resultArea = new JTextArea(10, 20);
        resultArea.setEditable(false);

        // Styling
        drawingPanel.setPreferredSize(new Dimension(280, 280));
        drawingPanel.setBorder(BorderFactory.createLineBorder(Color.GRAY));

        controlPanel.setLayout(new FlowLayout());
        controlPanel.add(clearButton);
        controlPanel.add(predictButton);

        add(drawingPanel, BorderLayout.CENTER);
        add(controlPanel, BorderLayout.SOUTH);
        add(new JScrollPane(resultArea), BorderLayout.EAST);

        // Event Listeners
        clearButton.addActionListener(e -> {
            drawingPanel.clear();
            resultArea.setText("");
        });

        predictButton.addActionListener(e -> predict());

        // Load Model
        loadModel();
    }

    private void loadModel() {
        try {
            network = MLOps.loadNN(MODEL_PATH);
            if (network == null) {
                resultArea.setText("Failed to load model from: " + MODEL_PATH + "\nCheck file path.");
            } else {
                resultArea.setText("Model Loaded Successfully!\nReady to draw.");
            }
        } catch (Exception e) {
            resultArea.setText("Error loading model: " + e.getMessage());
            e.printStackTrace();
        }
    }

    private void predict() {
        if (network == null) {
            resultArea.setText("Model not loaded.");
            return;
        }

        double[] input = drawingPanel.getPixelData();

        // Using ReLU as it is common for these networks (as seen in MLOps
        // heinitilizeWeights)
        // If the model was trained with Sigmoid, this might need changing.
        double[] output = MLOps.forwardPropagation(network, input, ActivationFunctions::reLU);

        displayResult(output);
    }

    private void displayResult(double[] output) {
        StringBuilder sb = new StringBuilder();
        int bestDigit = -1;
        double maxProb = -1;

        sb.append("Confidence Distribution:\n");
        sb.append("------------------------\n");

        for (int i = 0; i < output.length; i++) {
            sb.append(String.format("Digit %d: %.4f%%\n", i, output[i] * 100));
            if (output[i] > maxProb) {
                maxProb = output[i];
                bestDigit = i;
            }
        }

        sb.append("\nPrediction: " + bestDigit);
        sb.append(String.format("\nConfidence: %.2f%%", maxProb * 100));

        resultArea.setText(sb.toString());
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            new DrawApp().setVisible(true);
        });
    }

    // Inner class for the drawing canvas
    class DrawingPanel extends JPanel {
        private double[][] grid;
        private static final int GRID_SIZE = 28;
        private int lastX, lastY;

        public DrawingPanel() {
            setBackground(Color.BLACK);
            grid = new double[GRID_SIZE][GRID_SIZE];

            MouseAdapter mouseHandler = new MouseAdapter() {
                @Override
                public void mousePressed(MouseEvent e) {
                    lastX = e.getX();
                    lastY = e.getY();
                    draw(e.getX(), e.getY());
                }

                @Override
                public void mouseDragged(MouseEvent e) {
                    draw(e.getX(), e.getY());
                }
            };

            addMouseListener(mouseHandler);
            addMouseMotionListener(mouseHandler);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            int width = getWidth();
            int height = getHeight();
            double cellWidth = (double) width / GRID_SIZE;
            double cellHeight = (double) height / GRID_SIZE;

            for (int y = 0; y < GRID_SIZE; y++) {
                for (int x = 0; x < GRID_SIZE; x++) {
                    int val = (int) Math.min(255, Math.max(0, grid[y][x]));
                    g.setColor(new Color(val, val, val));
                    g.fillRect((int) (x * cellWidth), (int) (y * cellHeight),
                            (int) Math.ceil(cellWidth), (int) Math.ceil(cellHeight));
                }
            }

            // Optional: Draw grid lines for better visualization (can be removed if clean
            // look desired)
            g.setColor(new Color(30, 30, 30));
            for (int i = 0; i <= GRID_SIZE; i++) {
                g.drawLine((int) (i * cellWidth), 0, (int) (i * cellWidth), height);
                g.drawLine(0, (int) (i * cellHeight), width, (int) (i * cellHeight));
            }
        }

        private void draw(int x, int y) {
            // Bresenham-like interpolation or simple step interpolation
            double steps = Math.max(Math.abs(x - lastX), Math.abs(y - lastY));
            if (steps == 0) {
                addBrush(x, y);
            } else {
                double dx = (x - lastX) / steps;
                double dy = (y - lastY) / steps;
                for (int i = 0; i <= steps; i++) {
                    addBrush((int) (lastX + dx * i), (int) (lastY + dy * i));
                }
            }

            lastX = x;
            lastY = y;
            repaint();
        }

        private void addBrush(int px, int py) {
            double cellWidth = (double) getWidth() / GRID_SIZE;
            double cellHeight = (double) getHeight() / GRID_SIZE;

            // Convert pixel to grid coordinates
            double gx = px / cellWidth;
            double gy = py / cellHeight;

            // Brush radius in grid cells
            int radius = 1; // Increased radius to allow for smoother gradient

            for (int i = -radius; i <= radius; i++) {
                for (int j = -radius; j <= radius; j++) {
                    int cx = (int) gx + i;
                    int cy = (int) gy + j;

                    if (cx >= 0 && cx < GRID_SIZE && cy >= 0 && cy < GRID_SIZE) {
                        // Calculate distance from exact mouse position to cell center
                        double cellCenterX = cx + 0.5;
                        double cellCenterY = cy + 0.5;
                        double dist = Math.sqrt(Math.pow(gx - cellCenterX, 2) + Math.pow(gy - cellCenterY, 2));

                        // Gaussian-like falloff
                        double sigma = 1.0; // Controls the spread
                        double addedVal = 180 * Math.exp(-(dist * dist) / (2 * sigma * sigma));

                        // Threshold to avoid updating far pixels with negligible values
                        if (addedVal > 1) {
                            grid[cy][cx] += addedVal;
                            if (grid[cy][cx] > 255)
                                grid[cy][cx] = 255;
                        }
                    }
                }
            }
        }

        public void clear() {
            for (int y = 0; y < GRID_SIZE; y++) {
                for (int x = 0; x < GRID_SIZE; x++) {
                    grid[y][x] = 0;
                }
            }
            repaint();
        }

        public double[] getPixelData() {
            double[] pixels = new double[GRID_SIZE * GRID_SIZE];
            for (int y = 0; y < GRID_SIZE; y++) {
                for (int x = 0; x < GRID_SIZE; x++) {
                    pixels[y * GRID_SIZE + x] = Math.min(255, Math.max(0, grid[y][x])) / 255.0;
                }
            }
            return pixels;
        }
    }
}
