pub trait UnitLabel {
    fn symbol(&self) -> &'static str;
}

pub trait CoordinateLabel {
    fn quantity_name(&self) -> &'static str;
    fn quantity_symbol(&self) -> &'static str;
    fn unit_symbol(&self) -> &'static str;

    fn axis_label(&self) -> String {
        format!("{} ({})", self.quantity_symbol(), self.unit_symbol())
    }
}
