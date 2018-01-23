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
        class Vector2i
{

    static
            Vector2i Zero()
    {
        return new Vector2i(0, 0);
    }

    static
            Vector2i Unit()
    {
        return new Vector2i(1, 0);
    }

    static
            Vector2i One()
    {
        return new Vector2i(1, 1);
    }

    void Negate()
    {
        x *= -1;
        y *= -1;
    }

    void AddTo(Vector2i v2)
    {
        x += v2.x;
        y += v2.y;
    }

    Vector2i Add(Vector2i v2)
    {
        return new Vector2i(x + v2.x, y + v2.y);
    }

    void SubtractFrom(Vector2i v2)
    {
        x -= v2.x;
        y -= v2.y;
    }

    Vector2i Subtract(Vector2i v2)
    {
        return new Vector2i(x - v2.x, y - v2.y);
    }

    void AddTo(int i)
    {
        x += i;
        y += i;
    }

    Vector2i Add(int i)
    {
        return new Vector2i(x + i, y + i);
    }

    void SubtractFrom(int i)
    {
        x -= i;
        y -= i;
    }

    Vector2i Subtract(int i)
    {
        return new Vector2i(x - i, y - i);
    }

    void MultiplyBy(int d)
    {
        x *= d;
        y *= d;
    }

    Vector2i Multiply(int d)
    {
        return new Vector2i(x * d, y * d);
    }

    void DivideBy(int d)
    {
        x /= d;
        y /= d;
    }

    Vector2i Divide(int d)
    {
        return new Vector2i(x / d, y / d);
    }

    int Dot(Vector2i v2)
    {
        return (x * v2.x) + (y * v2.y);
    }

    int Cross(Vector2i v2)
    {
        return (x * v2.y) - (y * v2.x);
    }

    void HadamardBy(Vector2i v2)
    {
        x *= v2.x;
        y *= v2.y;
    }

    Vector2i Hadamard(Vector2i v2)
    {
        return new Vector2i(x * v2.x, y * v2.y);
    }

    void RotatateBy90By()
    {
        int tmp = y;
        y = x;
        x = -tmp;
    }

    Vector2i RotatateBy90()
    {
        return new Vector2i(-y, x);
    }

    void RotatateByMinus90By()
    {
        int tmp = x;
        x = y;
        y = -tmp;
    }

    Vector2i RotatateByMinus90()
    {
        return new Vector2i(y, -x);
    }

    Vector2i Copy()
    {
        return new Vector2i(x, y);
    }

    public
            int x, y;

    public
            Vector2i(int x, int y)
    {
        this.x = x;
        this.y = y;
    }

    public
            int LengthSq()
    {
        return (x * x) + (y * y);
    }

    public
            double Length()
    {
        return Math.sqrt(LengthSq());
    }

    public
            double DistanceSq(Vector2i v2)
    {
        v2 = Subtract(v2);
        return v2.LengthSq();
    }

    public
            double Distance(Vector2i v2)
    {
        v2 = Subtract(v2);
        return v2.Length();
    }

    Vector2 ToVector2()
    {
        return new Vector2(x, y);
    }

    @Override
    public
            String toString()
    {
        return "(" + x + ", " + y + ")";
    }
}
