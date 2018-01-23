/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package activecontourtracker;

/**
 *
 * @author edwardharry
 */
public
        class Vector2
{

    static
            Vector2 Zero()
    {
        return new Vector2(0, 0);
    }

    static
            Vector2 Unit()
    {
        return new Vector2(1, 0);
    }

    static
            Vector2 One()
    {
        return new Vector2(1, 1);
    }

    void Negate()
    {
        x *= -1;
        y *= -1;
    }

    void AddTo(Vector2 v2)
    {
        x += v2.x;
        y += v2.y;
    }

    Vector2 Add(Vector2 v2)
    {
        return new Vector2(x + v2.x, y + v2.y);
    }

    void SubtractFrom(Vector2 v2)
    {
        x -= v2.x;
        y -= v2.y;
    }

    Vector2 Subtract(Vector2 v2)
    {
        return new Vector2(x - v2.x, y - v2.y);
    }

    void AddTo(double d)
    {
        x += d;
        y += d;
    }

    Vector2 Add(double d)
    {
        return new Vector2(x + d, y + d);
    }

    void SubtractFrom(double d)
    {
        x -= d;
        y -= d;
    }

    Vector2 Subtract(double d)
    {
        return new Vector2(x - d, y - d);
    }

    void MultiplyBy(double d)
    {
        x *= d;
        y *= d;
    }

    Vector2 Multiply(double d)
    {
        return new Vector2(x * d, y * d);
    }

    void DivideBy(double d)
    {
        d = 1.0 / d;
        x *= d;
        y *= d;
    }

    Vector2 Divide(double d)
    {
        d = 1.0 / d;
        return new Vector2(x * d, y * d);
    }

    double Dot(Vector2 v2)
    {
        return (x * v2.x) + (y * v2.y);
    }

    double Cross(Vector2 v2)
    {
        return (x * v2.y) - (y * v2.x);
    }

    void HadamardBy(Vector2 v2)
    {
        x *= v2.x;
        y *= v2.y;
    }

    Vector2 Hadamard(Vector2 v2)
    {
        return new Vector2(x * v2.x, y * v2.y);
    }

    void RotatateBy90By()
    {
        double tmp = y;
        y = x;
        x = -tmp;
    }

    Vector2 RotatateBy90()
    {
        return new Vector2(-y, x);
    }

    void RotatateByMinus90By()
    {
        double tmp = x;
        x = y;
        y = -tmp;
    }

    Vector2 RotatateByMinus90()
    {
        return new Vector2(y, -x);
    }

    Vector2 Copy()
    {
        return new Vector2(x, y);
    }

    public
            double x, y;

    public
            Vector2(double x, double y)
    {
        this.x = x;
        this.y = y;
    }

    public
            double LengthSq()
    {
        return (x * x) + (y * y);
    }

    public
            double Length()
    {
        return Math.sqrt(LengthSq());
    }

    public
            double DistanceSq(Vector2 v2)
    {
        v2 = Subtract(v2);
        return v2.LengthSq();
    }

    public
            double Distance(Vector2 v2)
    {
        v2 = Subtract(v2);
        return v2.Length();
    }

    public
            void Normalise()
    {
        DivideBy(Length());
    }

    Vector2i ToVector2i()
    {
        return new Vector2i((int) (x + 0.5), (int) (y + 0.5));
    }

    Vector2i ToVector2i_floor()
    {
        return new Vector2i((int) x, (int) y);
    }

    @Override
    public
            String toString()
    {
        return "(" + x + ", " + y + ")";
    }
}
