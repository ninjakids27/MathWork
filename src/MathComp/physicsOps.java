public class physicsOps{
    private static final double G = 100;

  public static double velocity(double mass1, double radius){
    return (G*Math.sqrt(mass1/radius));
  }

  public static double period(double mass1, double radius){
    return (2*Math.PI*radius)/velocity(mass1, radius);
  }

  public static double centripitalAcceleration(double m1, double radius){
    return (Math.pow(velocity(m1, radius),2)/radius);
  }

  public static double gravitationalForce(double m1, double m2, double radius){
    return (m2*centripitalAcceleration(m1, radius));
  }

  public static void summary(double mass1, double mass2, double radius){
    System.out.println(ColorText.dataFormat("Velocity: " + velocity(mass1, radius)));
    System.out.println(ColorText.dataFormat("Period: " + period(mass1, radius)));
    System.out.println(ColorText.dataFormat("Centripital Acceleration: " + centripitalAcceleration(mass1, radius)));
    System.out.println(ColorText.dataFormat("Gravitational Force: " + gravitationalForce(mass1, mass2, radius)));
  }
}