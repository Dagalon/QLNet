namespace QLNet;

public class ShortRate1FYieldStructure<AffineModel1F> : YieldTermStructure
where AffineModel1F:OneFactorAffineModel
{
   #region Attributes

   private AffineModel1F model_;
   private double now_;
   private double r_t_;
   #endregion


   #region Constructors

   public ShortRate1FYieldStructure(Date refDate,  DayCounter dayCounter, Calendar calendar, AffineModel1F m)
   : base(refDate, calendar, dayCounter)
   {
      model_ = m;
      now_ = 0.0;
      double x0 = model_.dynamics().process().x0();
      r_t_ = model_.dynamics().shortRate(0.0, x0);
   }

   #endregion

   
   #region Methods

   public void updateState(double t, double r_t)
   {
      now_ = t;
      r_t_ = r_t;
   }

   #endregion


   public override Date maxDate()
   {
      return Date.maxDate();
   }

   protected override double discountImpl(double d)
   {
      return d > 0? model_.discountBond(now_, d, r_t_) : model_.discount(d);
   }
}
