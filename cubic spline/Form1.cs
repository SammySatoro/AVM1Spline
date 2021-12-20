using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using ZedGraph;


namespace cubic_spline
{
    

    public partial class Form1 : Form
    {

        double f(double x)     // вычисляемая функция
        {
            double tmp = Math.Sin(x);
             // double tmp = Math.Cos(x);

            return tmp;
        }

        public Form1()
        {
            InitializeComponent();
        }



        SplineTuple[] splines; // Сплайн
        SplineTuple1[] splins;

        // Структура, описывающая сплайн на каждом сегменте сетки
        private struct SplineTuple
        {
            public double a, b, c, d, x;
        }
        private struct SplineTuple1
        {
            public double a, b, c, d, x;
        }

        // Построение сплайна
        // x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
        // y - значения функции в узлах сетки
        // n - количество узлов сетки
        public void BuildSpline(double[] x, double[] y, int n)
        {
            // Инициализация массива сплайнов
            splines = new SplineTuple[n];
            for (int i = 0; i < n; ++i)
            {
                splines[i].x = x[i];
                splines[i].a = y[i];
            }
            splines[0].c = splines[n - 1].c = 0.0;

            // Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
            // Вычисление прогоночных коэффициентов - прямой ход метода прогонки
            double[] alpha = new double[n - 1];
            double[] beta = new double[n - 1];
            alpha[0] = beta[0] = 0.0;
            for (int i = 1; i < n - 1; ++i)
            {
                double hi = x[i] - x[i - 1];
                double hi1 = x[i + 1] - x[i];
                double A = hi;
                double C = 2.0 * (hi + hi1);
                double B = hi1;
                double F = 6.0 * ((y[i + 1] - y[i]) / hi1 - (y[i] - y[i - 1]) / hi);
                double z = (A * alpha[i - 1] + C);
                alpha[i] = -B / z;
                beta[i] = (F - A * beta[i - 1]) / z;
            }

            // Нахождение решения - обратный ход метода прогонки
            for (int i = n - 2; i > 0; --i)
            {
                splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
            }

            // По известным коэффициентам c[i] находим значения b[i] и d[i]
            for (int i = n - 1; i > 0; --i)
            {
                double hi = x[i] - x[i - 1];
                splines[i].d = (splines[i].c - splines[i - 1].c) / hi;
                splines[i].b = hi * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / hi;
            }
        }

        public void BuildSplineForAnother(double[] x, double[] y, int n)
        {
            // Инициализация массива сплайнов
            splins = new SplineTuple1[n];
            for (int i = 0; i < n; ++i)
            {
                splins[i].x = x[i];
                splins[i].a = y[i];
            }
            splins[0].c = splins[n - 1].c = 0.0;

            // Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
            // Вычисление прогоночных коэффициентов - прямой ход метода прогонки
            double[] alpha = new double[n - 1];
            double[] beta = new double[n - 1];
            alpha[0] = beta[0] = 0.0;
            for (int i = 1; i < n - 1; ++i)
            {
                double hi = x[i] - x[i - 1];
                double hi1 = x[i + 1] - x[i];
                double A = hi;
                double C = 2.0 * (hi + hi1);
                double B = hi1;
                double F = 6.0 * ((y[i + 1] - y[i]) / hi1 - (y[i] - y[i - 1]) / hi);
                double z = (A * alpha[i - 1] + C);
                alpha[i] = -B / z;
                beta[i] = (F - A * beta[i - 1]) / z;
            }

            // Нахождение решения - обратный ход метода прогонки
            for (int i = n - 2; i > 0; --i)
            {
                splins[i].c = alpha[i] * splins[i + 1].c + beta[i];
            }

            // По известным коэффициентам c[i] находим значения b[i] и d[i]
            for (int i = n - 1; i > 0; --i)
            {
                double hi = x[i] - x[i - 1];
                splins[i].d = (splins[i].c - splins[i - 1].c) / hi;
                splins[i].b = hi * (2.0 * splins[i].c + splins[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / hi;
            }
        }

        // Вычисление значения интерполированной функции в произвольной точке
        public double Interpolate(double x)
        {
            if (splines == null)
            {
                return double.NaN; // Если сплайны ещё не построены - возвращаем NaN
            }

            int n = splines.Length;
            SplineTuple s;

            if (x <= splines[0].x) // Если x меньше точки сетки x[0] - пользуемся первым эл-тов массива
            {
                s = splines[0];
            }
            else if (x >= splines[n - 1].x) // Если x больше точки сетки x[n - 1] - пользуемся последним эл-том массива
            {
                s = splines[n - 1];
            }
            else // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл-та массива
            {
                int i = 0;
                int j = n - 1;
                while (i + 1 < j)
                {
                    int k = i + (j - i) / 2;
                    if (x <= splines[k].x)
                    {
                        j = k;
                    }
                    else
                    {
                        i = k;
                    }
                }
                s = splines[j];
            }

            double dx = x - s.x;
            // Вычисляем значение сплайна в заданной точке по схеме Горнера (в принципе, "умный" компилятор применил бы схему Горнера сам, но ведь не все так умны, как кажутся)
            return s.a + (s.b + (s.c / 2.0 + s.d * dx / 6.0) * dx) * dx;
        }

        public double Interpolate2(double x)
        {
            if (splins == null)
            {
                return double.NaN; // Если сплайны ещё не построены - возвращаем NaN
            }

            int n = splins.Length;
            SplineTuple1 s;

            if (x <= splins[0].x) // Если x меньше точки сетки x[0] - пользуемся первым эл-тов массива
            {
                s = splins[0];
            }
            else if (x >= splins[n - 1].x) // Если x больше точки сетки x[n - 1] - пользуемся последним эл-том массива
            {
                s = splins[n - 1];
            }
            else // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл-та массива
            {
                int i = 0;
                int j = n - 1;
                while (i + 1 < j)
                {
                    int k = i + (j - i) / 2;
                    if (x <= splins[k].x)
                    {
                        j = k;
                    }
                    else
                    {
                        i = k;
                    }
                }
                s = splins[j];
            }

            double dx = x - s.x;
            // Вычисляем значение сплайна в заданной точке по схеме Горнера (в принципе, "умный" компилятор применил бы схему Горнера сам, но ведь не все так умны, как кажутся)
            return s.a + (s.b + (s.c / 2.0 + s.d * dx / 6.0) * dx) * dx;
        }

        public void Draw() // метод рисования для равносторонних узлов
        {
           double X = 2;
            int N = 13;
            double[] pointArrX = new double[N];
            double[] pointArrY = new double[N];


            GraphPane grap = zedGraphControl1.GraphPane; // графический объект

            grap.XAxis.Title.Text = "x";
            // Изменим текст по оси Y
            grap.YAxis.Title.Text = "y";
            // Изменим текст заголовка графика
            grap.Title.Text = "График функции y = sin(x)";

            grap.CurveList.Clear(); // очищаю график
            double Xmin = -10.0;
            double Xmax = 10.0;


            for (int i = 0; i < N; i++)
            { // создание равностоящих узлов
                pointArrX[i] = Xmin + (i * X);
                pointArrY[i] = f(Xmax + (i * X));
            }



            PointPairList point = new PointPairList();
            PointPairList basic = new PointPairList();
            double h = 1.0 / X;

            BuildSpline(pointArrX, pointArrY, N);

            for (double i = Xmin; i < Xmax; i += h)
            {
                point.Add(i, Interpolate(i)); // заполняю точки
            }


            


            for (double i = Xmin; i<Xmax; i+=h)
            {
              basic.Add(i, f(i));
            }

            LineItem line = grap.AddCurve("График сплайна для равностоящих", point, Color.Orange, SymbolType.Circle); // spline
           LineItem basic1 = grap.AddCurve("Начальная функция", basic, Color.DarkBlue, SymbolType.Diamond);
            zedGraphControl1.AxisChange();
            zedGraphControl1.Invalidate();
        }

        public void Draw2() // чебышевские узлы
        {
            double X = 2;
            int N = (int)X;
           
            double[] pointArrX1 = new double[N];
            double[] pointArrY1 = new double[N];

            GraphPane grap = zedGraphControl1.GraphPane; // графический объект

            grap.XAxis.Title.Text = "x";
            // Изменим текст по оси Y
            grap.YAxis.Title.Text = "y";
            // Изменим текст заголовка графика
            grap.Title.Text = "График функции y = sin(x)";

            grap.CurveList.Clear(); // очищаю график
            double Xmin = -10.0;
            double Xmax = 10.0;

            double pi = 3.14;
            for (int i = 0; i < N; i++)
            { // создание чебышевских узлов
                pointArrX1[i] = (Math.Cos(2 * i + 1) * pi) / (2 * N);
                pointArrY1[i] = f((Math.Cos(2 * i + 1) * pi) / (2 * N));
            }

            PointPairList anothernodes = new PointPairList();
            double h = 1.0 / X;

            BuildSplineForAnother(pointArrX1, pointArrY1, N);


            for (double i = Xmin; i < Xmax; i += h)
            {
                anothernodes.Add(i, Interpolate2(i)); // заполняю точки
            }

            zedGraphControl1.AxisChange();
            zedGraphControl1.Invalidate();
        }


        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void zedGraphControl1_Load(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
            Draw();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            Draw2();
        }
    }
}
