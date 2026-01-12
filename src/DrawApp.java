

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import MLComp.MLOps;
import MLComp.Neuron;
import MLComp.ActivationFunctions_Folder.ActivationFunctions;

public class DrawApp extends JFrame {

    private DrawingPanel drawingPanel;
    private JLabel predictionLabel;
    private JTextArea confidenceArea;
    private Neuron[][] network;
    private static final String MODEL_PATH = "Models//Optimus.ser";

    public DrawApp() {
        setTitle("MNIST Digit Recognizer - Neural Network Demo");
        setSize(800, 600);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout(10, 10));
        getContentPane().setBackground(new Color(30, 30, 40));

        // Main panel for drawing
        JPanel leftPanel = new JPanel(new BorderLayout(5, 5));
        leftPanel.setOpaque(false);

        JLabel drawLabel = new JLabel("Draw a digit (0-9)", SwingConstants.CENTER);
        drawLabel.setForeground(Color.WHITE);
        drawLabel.setFont(new Font("Arial", Font.BOLD, 16));
        leftPanel.add(drawLabel, BorderLayout.NORTH);

        drawingPanel = new DrawingPanel();
        drawingPanel.setPreferredSize(new Dimension(280, 280));
        drawingPanel.setBorder(BorderFactory.createLineBorder(new Color(100, 100, 150), 3));

        JPanel drawContainer = new JPanel(new FlowLayout(FlowLayout.CENTER));
        drawContainer.setOpaque(false);
        drawContainer.add(drawingPanel);
        leftPanel.add(drawContainer, BorderLayout.CENTER);

        // Control buttons
        JPanel controlPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 20, 10));
        controlPanel.setOpaque(false);

        JButton clearButton = createStyledButton("Clear", new Color(180, 60, 60));
        JButton predictButton = createStyledButton("Predict", new Color(60, 140, 60));

        controlPanel.add(clearButton);
        controlPanel.add(predictButton);
        leftPanel.add(controlPanel, BorderLayout.SOUTH);

        // Right panel for results
        JPanel rightPanel = new JPanel(new BorderLayout(10, 10));
        rightPanel.setOpaque(false);
        rightPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 20));

        // Large prediction display
        predictionLabel = new JLabel("?", SwingConstants.CENTER);
        predictionLabel.setFont(new Font("Arial", Font.BOLD, 150));
        predictionLabel.setForeground(new Color(100, 200, 100));
        predictionLabel.setPreferredSize(new Dimension(200, 200));
        predictionLabel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(new Color(100, 100, 150), 2),
                BorderFactory.createEmptyBorder(20, 20, 20, 20)));

        JPanel predictionPanel = new JPanel(new BorderLayout());
        predictionPanel.setOpaque(false);
        JLabel predTitle = new JLabel("Prediction", SwingConstants.CENTER);
        predTitle.setForeground(Color.WHITE);
        predTitle.setFont(new Font("Arial", Font.BOLD, 18));
        predictionPanel.add(predTitle, BorderLayout.NORTH);
        predictionPanel.add(predictionLabel, BorderLayout.CENTER);

        rightPanel.add(predictionPanel, BorderLayout.NORTH);

        // Confidence distribution
        JLabel confLabel = new JLabel("Confidence Distribution", SwingConstants.CENTER);
        confLabel.setForeground(Color.WHITE);
        confLabel.setFont(new Font("Arial", Font.BOLD, 14));

        confidenceArea = new JTextArea(12, 25);
        confidenceArea.setEditable(false);
        confidenceArea.setFont(new Font("Monospaced", Font.PLAIN, 14));
        confidenceArea.setBackground(new Color(40, 40, 50));
        confidenceArea.setForeground(new Color(200, 200, 200));
        confidenceArea.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        JPanel confPanel = new JPanel(new BorderLayout(5, 5));
        confPanel.setOpaque(false);
        confPanel.add(confLabel, BorderLayout.NORTH);
        confPanel.add(new JScrollPane(confidenceArea), BorderLayout.CENTER);

        rightPanel.add(confPanel, BorderLayout.CENTER);

        // Status label
        JLabel statusLabel = new JLabel("Model: " + MODEL_PATH, SwingConstants.CENTER);
        statusLabel.setForeground(new Color(150, 150, 150));
        statusLabel.setFont(new Font("Arial", Font.ITALIC, 11));
        rightPanel.add(statusLabel, BorderLayout.SOUTH);

        add(leftPanel, BorderLayout.CENTER);
        add(rightPanel, BorderLayout.EAST);

        // Event Listeners
        clearButton.addActionListener(e -> {
            drawingPanel.clear();
            predictionLabel.setText("?");
            predictionLabel.setForeground(new Color(100, 200, 100));
            confidenceArea.setText("");
        });

        predictButton.addActionListener(e -> predict());

        // Load Model
        loadModel();
    }

    private JButton createStyledButton(String text, Color bgColor) {
        JButton button = new JButton(text);
        button.setFont(new Font("Arial", Font.BOLD, 14));
        button.setBackground(bgColor);
        button.setForeground(Color.WHITE);
        button.setFocusPainted(false);
        button.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createRaisedBevelBorder(),
                BorderFactory.createEmptyBorder(8, 25, 8, 25)));
        button.setCursor(new Cursor(Cursor.HAND_CURSOR));
        return button;
    }

    private void loadModel() {
        try {
            network = MLOps.loadNN(MODEL_PATH);
            if (network == null) {
                confidenceArea.setText("Failed to load model!\nCheck: " + MODEL_PATH);
                predictionLabel.setText("!");
                predictionLabel.setForeground(Color.RED);
            } else {
                confidenceArea.setText("Model loaded successfully!\n\nDraw a digit and click Predict.");
            }
        } catch (Exception e) {
            confidenceArea.setText("Error: " + e.getMessage());
            predictionLabel.setText("!");
            predictionLabel.setForeground(Color.RED);
            e.printStackTrace();
        }
    }

    private void predict() {
        if (network == null) {
            confidenceArea.setText("Model not loaded.");
            return;
        }

        double[] input = drawingPanel.getPixelData();
        double[] output = MLOps.forwardPropagation(network, input, ActivationFunctions::reLU);
        displayResult(output);
    }

    private void displayResult(double[] output) {
        StringBuilder sb = new StringBuilder();
        int bestDigit = -1;
        double maxProb = -1;

        for (int i = 0; i < output.length; i++) {
            if (output[i] > maxProb) {
                maxProb = output[i];
                bestDigit = i;
            }
        }

        // Update large prediction display
        predictionLabel.setText(String.valueOf(bestDigit));
        if (maxProb > 0.8) {
            predictionLabel.setForeground(new Color(100, 200, 100)); // Green for confident
        } else if (maxProb > 0.5) {
            predictionLabel.setForeground(new Color(200, 200, 100)); // Yellow for moderate
        } else {
            predictionLabel.setForeground(new Color(200, 100, 100)); // Red for uncertain
        }

        // Build confidence bars
        for (int i = 0; i < output.length; i++) {
            double pct = output[i] * 100;
            int barLength = (int) (pct / 5); // Scale to fit
            StringBuilder barBuilder = new StringBuilder();
            for (int b = 0; b < Math.max(0, barLength); b++) {
                barBuilder.append("█");
            }
            String bar = barBuilder.toString();
            String marker = (i == bestDigit) ? " ◄" : "";
            sb.append(String.format("%d: %6.2f%% %s%s\n", i, pct, bar, marker));
        }

        sb.append(String.format("\nConfidence: %.1f%%", maxProb * 100));

        confidenceArea.setText(sb.toString());
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (Exception ignored) {
            }
            new DrawApp().setVisible(true);
        });
    }

    // Inner class for the drawing canvas
    // Uses higher resolution grid (CANVAS_SIZE) for smoother drawing, then
    // downscales to OUTPUT_SIZE for NN
    class DrawingPanel extends JPanel {
        private double[][] grid;
        private static final int CANVAS_SIZE = 112; // High-res canvas (4x upscale)
        private static final int OUTPUT_SIZE = 28; // Final size for neural network
        private static final int SCALE_FACTOR = CANVAS_SIZE / OUTPUT_SIZE; // 4x
        private int lastX = -1, lastY = -1;

        public DrawingPanel() {
            setBackground(Color.BLACK);
            grid = new double[CANVAS_SIZE][CANVAS_SIZE];

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

                @Override
                public void mouseReleased(MouseEvent e) {
                    lastX = -1;
                    lastY = -1;
                }
            };

            addMouseListener(mouseHandler);
            addMouseMotionListener(mouseHandler);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);

            int width = getWidth();
            int height = getHeight();
            double cellWidth = (double) width / CANVAS_SIZE;
            double cellHeight = (double) height / CANVAS_SIZE;

            for (int y = 0; y < CANVAS_SIZE; y++) {
                for (int x = 0; x < CANVAS_SIZE; x++) {
                    int val = (int) Math.min(255, Math.max(0, grid[y][x]));
                    g.setColor(new Color(val, val, val));
                    g.fillRect((int) (x * cellWidth), (int) (y * cellHeight),
                            (int) Math.ceil(cellWidth) + 1, (int) Math.ceil(cellHeight) + 1);
                }
            }
        }

        private void draw(int x, int y) {
            if (lastX < 0 || lastY < 0) {
                addBrush(x, y);
            } else {
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
            }

            lastX = x;
            lastY = y;
            repaint();
        }

        private void addBrush(int px, int py) {
            double cellWidth = (double) getWidth() / CANVAS_SIZE;
            double cellHeight = (double) getHeight() / CANVAS_SIZE;

            double gx = px / cellWidth;
            double gy = py / cellHeight;

            // Larger brush radius for smoother strokes (scaled to higher res canvas)
            int radius = 2;
            double sigma = 3; // Gaussian spread for soft edges

            for (int i = -radius; i <= radius; i++) {
                for (int j = -radius; j <= radius; j++) {
                    int cx = (int) gx + i;
                    int cy = (int) gy + j;

                    if (cx >= 0 && cx < CANVAS_SIZE && cy >= 0 && cy < CANVAS_SIZE) {
                        double cellCenterX = cx + 0.5;
                        double cellCenterY = cy + 0.5;
                        double dist = Math.sqrt(Math.pow(gx - cellCenterX, 2) + Math.pow(gy - cellCenterY, 2));

                        // Gaussian falloff for smooth gradient edges
                        double addedVal = 200 * Math.exp(-(dist * dist) / (2 * sigma * sigma));

                        if (addedVal > 0.5) {
                            grid[cy][cx] += addedVal;
                            if (grid[cy][cx] > 255) {
                                grid[cy][cx] = 255;
                            }
                        }
                    }
                }
            }
        }

        public void clear() {
            for (int y = 0; y < CANVAS_SIZE; y++) {
                for (int x = 0; x < CANVAS_SIZE; x++) {
                    grid[y][x] = 0;
                }
            }
            lastX = -1;
            lastY = -1;
            repaint();
        }

        /**
         * Downscales the high-res canvas to 28x28 using bilinear interpolation
         * and returns normalized pixel data for the neural network
         */
        public double[] getPixelData() {
            double[] pixels = new double[OUTPUT_SIZE * OUTPUT_SIZE];

            // Bilinear downsampling from CANVAS_SIZE to OUTPUT_SIZE
            for (int y = 0; y < OUTPUT_SIZE; y++) {
                for (int x = 0; x < OUTPUT_SIZE; x++) {
                    double sum = 0;
                    // Average the SCALE_FACTOR x SCALE_FACTOR block
                    for (int dy = 0; dy < SCALE_FACTOR; dy++) {
                        for (int dx = 0; dx < SCALE_FACTOR; dx++) {
                            int srcY = y * SCALE_FACTOR + dy;
                            int srcX = x * SCALE_FACTOR + dx;
                            sum += grid[srcY][srcX];
                        }
                    }
                    double avgVal = sum / (SCALE_FACTOR * SCALE_FACTOR);
                    pixels[y * OUTPUT_SIZE + x] = Math.min(255, Math.max(0, avgVal)) / 255.0;
                }
            }
            return pixels;
        }
    }
}
